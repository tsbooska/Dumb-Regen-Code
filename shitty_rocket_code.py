"""
Liquid rocket engine preliminary design tool.

Predicts delivered performance and generates chamber + nozzle geometry from:
    chamber pressure, thrust, propellants, O/F ratio, L*, and efficiency factors.

Thermochemistry comes from NASA CEA via RocketCEA (pip install rocketcea).
Everything else (sizing, losses, contour) is computed here.

Method
------
1. CEA (equilibrium or frozen) at Pc, O/F, and the expansion ratio that gives the
   design exit pressure  ->  ideal c*, Cf_vac, Tc, gamma, MW, Pe, Me, ...
2. Losses:
     c*_del   = eta_cstar * c*_ideal                         (injector / combustion)
     Cf_vac   = eta_div * eta_nozzle * Cf_vac_ideal          (divergence, BL, kinetics)
     eta_div  = (1 + cos(alpha)) / 2                         (conical divergence factor)
     Cf(Pa)   = Cf_vac - eps * Pa / Pc
     Isp      = c*_del * Cf / g0
3. Sizing:   At = F / (Pc * Cf(Pa)),   mdot = Pc * At / c*_del
4. Chamber:  Vc = L* * At  (volume from injector face to throat, incl. convergent)
             Contraction ratio from input or Humble's correlation
             Ac/At = 8.0 * Dt[cm]^-0.6 + 1.25
             Cylinder + straight convergent cone (frustum), sharp corners, no fillets
5. Nozzle:   Straight conical divergent section from a sharp throat.

Coordinate system: x = 0 at the throat, negative upstream (chamber), positive downstream.

Usage examples
--------------
    python shitty_rocket_code.py                                     # defaults: LOX/RP-1, 5 kN, 20 bar
    python shitty_rocket_code.py --ox LOX --fuel CH4 --pc 50 --thrust 25000 --of 3.2 --lstar 1.0
    python shitty_rocket_code.py --ox N2O --fuel Ethanol --of-sweep 2 6 0.25 --optimize-of
    python shitty_rocket_code.py --pamb 0 --eps 60 --cone-angle 12 --plot

Common RocketCEA propellant names:
    oxidizers: LOX, GOX, N2O4, N2O, H2O2, IRFNA, MON25, F2
    fuels:     RP-1, LH2, GH2, CH4, LCH4_NASA, Ethanol, Methanol, MMH, UDMH, N2H4, A50, JetA, Propane
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field

import numpy as np

try:
    from rocketcea.cea_obj_w_units import CEA_Obj
except ImportError:  # pragma: no cover
    sys.exit("RocketCEA is required:  pip install rocketcea")

G0 = 9.80665          # m/s^2
BAR = 1.0e5           # Pa
R_UNIV = 8314.462618  # J/(kmol K)


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------
@dataclass
class EngineInputs:
    # Operating point
    pc_bar: float = 20.0                 # chamber pressure (injector-end stagnation), bar
    thrust_n: float = 1800             # thrust at the design ambient pressure, N
    pamb_bar: float = 1.01325            # design ambient pressure, bar (0 = vacuum)
    pe_bar: float | None = None          # design exit pressure, bar (default = pamb -> optimum)
    eps: float | None = None             # expansion ratio override (takes priority over pe)

    # Propellants
    ox: str = "LOX"
    fuel: str = "RP-1"
    of_ratio: float = 1.3                # CORE O/F (core oxidizer / core fuel); film fuel not included
    film_pct_core: float = 0.0           # fuel film flow as % of the CORE mass flow (ox + core fuel); 0 = none
    of_sweep: tuple[float, float, float] | None = None   # (min, max, step)
    optimize_of: bool = False            # use the O/F with max delivered Isp from the sweep
    frozen: bool = False                 # frozen composition (from chamber) instead of equilibrium
    frozen_at_throat: bool = False       # freeze at the throat instead of the chamber

    # Chamber
    lstar_m: float = 1.1                 # characteristic length L*, m
    contraction_ratio: float | None = 4   # Ac/At; None -> Humble correlation
    conv_half_angle_deg: float = 30.0    # convergent cone half angle

    # Nozzle
    cone_half_angle_deg: float = 15.0    # divergent cone half angle

    # Efficiencies
    eta_cstar: float = 0.95              # injector / combustion efficiency (c* efficiency)
    eta_nozzle: float = 0.98             # nozzle losses other than divergence (BL, kinetics)
    eta_divergence: float | None = None  # override computed divergence efficiency

    # Output
    csv_path: str | None = "engine_contour.csv"
    plot: bool = False
    plot_path: str | None = "engine_contour.png"


# --------------------------------------------------------------------------------------
# Thermochemistry / performance
# --------------------------------------------------------------------------------------
def make_cea(ox: str, fuel: str) -> CEA_Obj:
    return CEA_Obj(
        oxName=ox, fuelName=fuel,
        pressure_units="bar", temperature_units="K", cstar_units="m/s",
        isp_units="sec", sonic_velocity_units="m/s", enthalpy_units="J/kg",
        density_units="kg/m^3", specific_heat_units="J/kg-K",
    )


def evaluate_performance(inp: EngineInputs, cea: CEA_Obj, of: float) -> dict:
    """Ideal (CEA) and delivered performance at one O/F. Does not size the engine."""
    pc = inp.pc_bar
    fz = dict(frozen=int(inp.frozen or inp.frozen_at_throat), frozenAtThroat=int(inp.frozen_at_throat))

    # Expansion ratio
    if inp.eps is not None:
        eps = inp.eps
    else:
        pe = inp.pe_bar if inp.pe_bar is not None else inp.pamb_bar
        if pe <= 0:
            raise ValueError("Vacuum design: specify --eps or a non-zero --pe.")
        if pe >= pc:
            raise ValueError("Exit pressure must be below chamber pressure.")
        eps = float(cea.get_eps_at_PcOvPe(Pc=pc, MR=of, PcOvPe=pc / pe, **fz))

    isp_vac_ideal, cstar_ideal, tc = (float(v) for v in cea.get_IvacCstrTc(Pc=pc, MR=of, eps=eps, **fz))
    cf_vac_ideal = isp_vac_ideal * G0 / cstar_ideal
    pe_bar = pc / float(cea.get_PcOvPe(Pc=pc, MR=of, eps=eps, **fz))
    mw_c, gam_c = (float(v) for v in cea.get_Chamber_MolWt_gamma(Pc=pc, MR=of, eps=eps))
    mw_t, gam_t = (float(v) for v in cea.get_Throat_MolWt_gamma(Pc=pc, MR=of, eps=eps, frozen=fz["frozen"]))
    temps = [float(t) for t in cea.get_Temperatures(Pc=pc, MR=of, eps=eps, **fz)]
    rhos = [float(r) for r in cea.get_Densities(Pc=pc, MR=of, eps=eps, **fz)]
    mach_e = float(cea.get_MachNumber(Pc=pc, MR=of, eps=eps, **fz))

    # Losses
    if inp.eta_divergence is not None:
        eta_div = inp.eta_divergence
    else:
        # Conical source-flow divergence factor
        eta_div = 0.5 * (1.0 + math.cos(math.radians(inp.cone_half_angle_deg)))
    eta_cf = eta_div * inp.eta_nozzle

    cstar = inp.eta_cstar * cstar_ideal
    cf_vac = eta_cf * cf_vac_ideal
    pa_pc = inp.pamb_bar / pc
    cf = cf_vac - eps * pa_pc
    cf_ideal = cf_vac_ideal - eps * pa_pc
    cf_sl = cf_vac - eps * 1.01325 / pc

    return dict(
        of=of, eps=eps, pe_bar=pe_bar, mach_e=mach_e,
        tc=tc, t_throat=temps[1], t_exit=temps[2],
        rho_c=rhos[0], mw_c=mw_c, gam_c=gam_c, mw_t=mw_t, gam_t=gam_t,
        cstar_ideal=cstar_ideal, cf_vac_ideal=cf_vac_ideal, cf_ideal=cf_ideal,
        isp_vac_ideal=isp_vac_ideal, isp_ideal=cstar_ideal * cf_ideal / G0,
        eta_div=eta_div, eta_cf=eta_cf,
        cstar=cstar, cf_vac=cf_vac, cf=cf, cf_sl=cf_sl,
        isp=cstar * cf / G0, isp_vac=cstar * cf_vac / G0, isp_sl=cstar * cf_sl / G0,
    )


def sweep_of(inp: EngineInputs, cea: CEA_Obj) -> list[dict]:
    lo, hi, step = inp.of_sweep
    results = []
    for of in np.arange(lo, hi + 0.5 * step, step):
        try:
            results.append(evaluate_performance(inp, cea, float(of)))
        except Exception as exc:  # CEA can fail at extreme mixture ratios
            print(f"  O/F {of:.3f}: CEA failed ({exc})")
    return results


# --------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------
def humble_contraction_ratio(dt_m: float) -> float:
    """Ac/At = 8.0 * Dt^-0.6 + 1.25 with Dt in cm (Humble, Space Propulsion Analysis & Design)."""
    return 8.0 * (dt_m * 100.0) ** -0.6 + 1.25


def chamber_contour(inp: EngineInputs, rt: float, rc: float, vc: float, warnings: list[str]) -> dict:
    """Injector face -> throat: cylinder + straight convergent cone. Returns contour and dimensions (SI)."""
    th = math.radians(inp.conv_half_angle_deg)

    # Convergent frustum from Rc down to a sharp throat at x = 0
    l_conv = (rc - rt) / math.tan(th)
    v_conv = math.pi * l_conv / 3.0 * (rc ** 2 + rc * rt + rt ** 2)

    l_cyl = (vc - v_conv) / (math.pi * rc ** 2)
    if l_cyl < 0:
        warnings.append("L* volume is smaller than the convergent section alone: cylindrical length set "
                        "to 0. Increase L* or reduce the contraction ratio / increase the convergent angle.")
        l_cyl = 0.0
    x_inj = -l_conv - l_cyl

    x = np.array([x_inj, -l_conv, 0.0])
    r = np.array([rc, rc, rt])
    return dict(x=x, r=r, x_inj=x_inj, l_cyl=l_cyl, l_conv=l_conv, v_conv=v_conv,
                v_cyl=math.pi * rc ** 2 * l_cyl)


def nozzle_contour(inp: EngineInputs, rt: float, re: float) -> dict:
    """Throat -> exit: straight conical divergent section."""
    ln = (re - rt) / math.tan(math.radians(inp.cone_half_angle_deg))
    return dict(x=np.array([0.0, ln]), r=np.array([rt, re]), length=ln)


# --------------------------------------------------------------------------------------
# Top-level design
# --------------------------------------------------------------------------------------
@dataclass
class EngineDesign:
    inputs: EngineInputs
    perf: dict
    sizing: dict
    chamber: dict
    nozzle: dict
    sweep: list[dict] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def contour(self) -> tuple[np.ndarray, np.ndarray]:
        return (np.concatenate([self.chamber["x"], self.nozzle["x"][1:]]),
                np.concatenate([self.chamber["r"], self.nozzle["r"][1:]]))


def design_engine(inp: EngineInputs) -> EngineDesign:
    warnings: list[str] = []
    cea = make_cea(inp.ox, inp.fuel)

    sweep = sweep_of(inp, cea) if inp.of_sweep else []
    of = inp.of_ratio
    if inp.optimize_of:
        if not sweep:
            raise ValueError("--optimize-of requires --of-sweep.")
        of = max(sweep, key=lambda p: p["isp"])["of"]

    perf = evaluate_performance(inp, cea, of)
    if perf["cf"] <= 0:
        raise ValueError("Thrust coefficient at the design ambient pressure is <= 0; "
                         "nozzle is far too overexpanded for this Pc.")

    # Sizing
    pc = inp.pc_bar * BAR
    at = inp.thrust_n / (pc * perf["cf"])
    rt = math.sqrt(at / math.pi)
    ae = at * perf["eps"]
    re = math.sqrt(ae / math.pi)
    mdot = pc * at / perf["cstar"]
    mdot_ox = mdot * of / (1 + of)
    mdot_f = mdot / (1 + of)

    cr = inp.contraction_ratio or humble_contraction_ratio(2 * rt)
    ac = cr * at
    rc = math.sqrt(ac / math.pi)
    vc = inp.lstar_m * at

    chamber = chamber_contour(inp, rt, rc, vc, warnings)
    nozzle = nozzle_contour(inp, rt, re)

    # Chamber gas residence time  t_s = m_gas / mdot = rho_c * Vc / mdot
    t_res = perf["rho_c"] * vc / mdot

    # Mach number at the chamber/convergent junction (isentropic, frozen gamma), indicates
    # Rayleigh (heat-addition) stagnation-pressure loss when > ~0.3
    g = perf["gam_c"]
    mc = _subsonic_mach(cr, g)

    sizing = dict(
        at=at, rt=rt, dt=2 * rt, ae=ae, re=re, de=2 * re, ac=ac, rc=rc, dc=2 * rc, cr=cr, vc=vc,
        mdot=mdot, mdot_ox=mdot_ox, mdot_f=mdot_f, t_res=t_res, mach_c=mc,
        thrust_vac=perf["cf_vac"] * pc * at, thrust_sl=perf["cf_sl"] * pc * at,
        length_total=nozzle["length"] - chamber["x_inj"],
    )

    # Fuel film: a % of the CORE flow (non-standard convention). The core (mdot, mdot_ox, mdot_f, O/F)
    # is sized as without film; the film adds mass flow but no thrust (throat blockage neglected).
    mdot_film = inp.film_pct_core / 100.0 * mdot
    sizing.update(mdot_core=mdot, mdot_film=mdot_film, mdot_total=mdot + mdot_film,
                  mdot_f_total=mdot_f + mdot_film)
    perf["of_overall"] = mdot_ox / (mdot_f + mdot_film)
    for key, thrust in (("isp_film", inp.thrust_n), ("isp_vac_film", sizing["thrust_vac"]),
                        ("isp_sl_film", sizing["thrust_sl"])):
        perf[key] = thrust / ((mdot + mdot_film) * G0)
    if mdot_film > 0:
        warnings.append("Film cooling: the film is assumed to add mass flow but no thrust, and its blockage of "
                        "the throat is neglected (throat sized on the core flow).")

    # Sanity checks
    pe, pa = perf["pe_bar"], inp.pamb_bar
    if pa > 0 and pe < 0.4 * pa:
        warnings.append(f"Pe = {pe:.3f} bar < 0.4 Pa: likely flow separation at the design point "
                        "(Summerfield criterion).")
    if pe * BAR < 0.4 * 1.01325 * BAR and pa != 1.01325:
        warnings.append(f"Nozzle would likely separate at sea level (Pe = {pe:.3f} bar).")
    if mc > 0.3:
        warnings.append(f"Chamber Mach {mc:.2f} > 0.3: significant Rayleigh Pc loss; "
                        "consider a larger contraction ratio.")
    if t_res < 1e-3:
        warnings.append(f"Residence time {t_res * 1e3:.2f} ms is short; check L* for these propellants.")

    return EngineDesign(inp, perf, sizing, chamber, nozzle, sweep, warnings)


def _subsonic_mach(area_ratio: float, g: float) -> float:
    """Solve A/A* = f(M) on the subsonic branch by bisection."""
    def a_ratio(m):
        return (1 / m) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * m * m)) ** ((g + 1) / (2 * (g - 1)))
    lo, hi = 1e-6, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if a_ratio(mid) > area_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _supersonic_mach(area_ratio: float, g: float) -> float:
    """Solve A/A* = f(M) on the supersonic branch by bisection."""
    def a_ratio(m):
        return (1 / m) * ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * m * m)) ** ((g + 1) / (2 * (g - 1)))
    lo, hi = 1.0, 50.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if a_ratio(mid) < area_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def print_report(d: EngineDesign) -> None:
    i, p, s, c, n = d.inputs, d.perf, d.sizing, d.chamber, d.nozzle
    mm = 1e3
    mode = "frozen at throat" if i.frozen_at_throat else ("frozen" if i.frozen else "shifting equilibrium")

    if d.sweep:
        print("\nO/F SWEEP (delivered values at the design ambient pressure)")
        print(f"  {'O/F':>6} {'eps':>7} {'Tc [K]':>8} {'c* [m/s]':>9} {'Cf':>7} {'Isp [s]':>8} {'Isp_vac':>8}")
        best = max(d.sweep, key=lambda q: q["isp"])
        for q in d.sweep:
            mark = "  <- max Isp" if q is best else ""
            print(f"  {q['of']:6.2f} {q['eps']:7.2f} {q['tc']:8.0f} {q['cstar']:9.1f} "
                  f"{q['cf']:7.4f} {q['isp']:8.1f} {q['isp_vac']:8.1f}{mark}")

    print("\n" + "=" * 70)
    print(f" ENGINE DESIGN: {i.ox} / {i.fuel}   ({mode})")
    print("=" * 70)
    print(f" Chamber pressure          {i.pc_bar:10.2f} bar")
    print(f" Design thrust             {i.thrust_n:10.1f} N   at Pa = {i.pamb_bar:.4f} bar")
    print(f" O/F ratio                 {p['of']:10.3f}")
    print(f" L*                        {i.lstar_m:10.3f} m")
    print(f" eta c* / eta nozzle / eta div  {i.eta_cstar:.3f} / {i.eta_nozzle:.3f} / {p['eta_div']:.4f}")

    print("\n THERMOCHEMISTRY (CEA)")
    print(f"  Tc / Tt / Te             {p['tc']:8.0f} / {p['t_throat']:.0f} / {p['t_exit']:.0f} K")
    print(f"  MW, gamma (chamber)      {p['mw_c']:8.3f} kg/kmol, {p['gam_c']:.4f}")
    print(f"  MW, gamma (throat)       {p['mw_t']:8.3f} kg/kmol, {p['gam_t']:.4f}")
    print(f"  Exit pressure            {p['pe_bar']:8.4f} bar")
    print(f"  Exit Mach                {p['mach_e']:8.3f}")

    print("\n PERFORMANCE                ideal      delivered")
    print(f"  c* [m/s]               {p['cstar_ideal']:9.1f}   {p['cstar']:9.1f}")
    print(f"  Cf (design Pa)         {p['cf_ideal']:9.4f}   {p['cf']:9.4f}")
    print(f"  Cf (vacuum)            {p['cf_vac_ideal']:9.4f}   {p['cf_vac']:9.4f}")
    print(f"  Isp (design Pa) [s]    {p['isp_ideal']:9.1f}   {p['isp']:9.1f}")
    print(f"  Isp (vacuum) [s]       {p['isp_vac_ideal']:9.1f}   {p['isp_vac']:9.1f}")
    print(f"  Isp (sea level) [s]                {p['isp_sl']:9.1f}")
    print(f"  Thrust vac / SL [N]                {s['thrust_vac']:9.1f} / {s['thrust_sl']:.1f}")

    print("\n MASS FLOW")
    print(f"  Total                    {s['mdot']:10.4f} kg/s")
    print(f"  Oxidizer                 {s['mdot_ox']:10.4f} kg/s")
    print(f"  Fuel                     {s['mdot_f']:10.4f} kg/s")
    if s["mdot_film"] > 0:
        print("  (above: core flow; O/F is the core O/F)")
        print(f"  Film fuel                {s['mdot_film']:10.4f} kg/s   ({i.film_pct_core:g} % of core flow)")
        print(f"  Total fuel / total flow  {s['mdot_f_total']:10.4f} / {s['mdot_total']:.4f} kg/s   "
              f"(overall O/F {p['of_overall']:.3f})")
        print(f"  Isp with film (design Pa){p['isp_film']:10.1f} s   (vacuum {p['isp_vac_film']:.1f} s; "
              "film adds mass, no thrust)")
    print(f"  Residence time           {s['t_res'] * 1e3:10.3f} ms")

    print("\n CHAMBER GEOMETRY")
    print(f"  Chamber diameter         {s['dc'] * mm:10.2f} mm")
    print(f"  Contraction ratio Ac/At  {s['cr']:10.3f}{'  (Humble correlation)' if not i.contraction_ratio else ''}")
    print(f"  Chamber Mach (inlet)     {s['mach_c']:10.4f}")
    print(f"  Chamber volume (L*·At)   {s['vc'] * 1e6:10.2f} cm^3")
    print(f"    cylindrical part       {c['v_cyl'] * 1e6:10.2f} cm^3")
    print(f"    convergent part        {c['v_conv'] * 1e6:10.2f} cm^3")
    print(f"  Cylindrical length       {c['l_cyl'] * mm:10.2f} mm")
    print(f"  Convergent length        {c['l_conv'] * mm:10.2f} mm")
    print(f"  Convergent half angle    {i.conv_half_angle_deg:10.2f} deg")

    print("\n THROAT & NOZZLE GEOMETRY")
    print(f"  Throat diameter          {s['dt'] * mm:10.2f} mm   (At = {s['at'] * 1e6:.2f} mm^2)")
    print(f"  Exit diameter            {s['de'] * mm:10.2f} mm   (Ae = {s['ae'] * 1e6:.1f} mm^2)")
    print(f"  Expansion ratio Ae/At    {p['eps']:10.3f}")
    print(f"  Divergent half angle     {i.cone_half_angle_deg:10.2f} deg")
    print(f"  Nozzle length            {n['length'] * mm:10.2f} mm")
    print(f"\n  TOTAL LENGTH (injector face -> exit)  {s['length_total'] * mm:.2f} mm")

    if d.warnings:
        print("\n WARNINGS")
        for w in d.warnings:
            print(f"  ! {w}")
    print()


def write_csv(d: EngineDesign, path: str) -> None:
    x, r = d.contour
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x_mm", "r_mm"])
        for xi, ri in zip(x, r):
            w.writerow([f"{xi * 1e3:.4f}", f"{ri * 1e3:.4f}"])
    print(f"Contour written to {path}  ({len(x)} points, x = 0 at throat)")


def plot_contour(d: EngineDesign, path: str | None, show: bool) -> None:
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, r = d.contour
    x, r = x * 1e3, r * 1e3
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(x, r, "k", lw=1.8)
    ax.plot(x, -r, "k", lw=1.8)
    ax.fill_between(x, r, -r, color="tab:orange", alpha=0.12)
    ax.axhline(0, color="gray", lw=0.6, ls="-.")
    ax.axvline(0, color="tab:red", lw=0.8, ls="--", label="throat")
    ax.axvline(d.chamber["x_inj"] * 1e3, color="tab:blue", lw=0.8, ls="--", label="injector face")
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    p, s, i = d.perf, d.sizing, d.inputs
    ax.set_title(f"{i.ox}/{i.fuel}  Pc={i.pc_bar:g} bar  F={i.thrust_n:g} N  O/F={p['of']:.2f}  "
                 f"eps={p['eps']:.1f}  Isp={p['isp']:.0f} s  Dt={s['dt'] * 1e3:.1f} mm")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150)
        print(f"Plot written to {path}")
    if show:
        plt.show()


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def parse_args(argv=None) -> tuple[EngineInputs, bool]:
    d = EngineInputs()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a = ap.add_argument
    a("--pc", type=float, default=d.pc_bar, help="chamber pressure [bar]")
    a("--thrust", type=float, default=d.thrust_n, help="thrust at design ambient pressure [N]")
    a("--pamb", type=float, default=d.pamb_bar, help="design ambient pressure [bar], 0 = vacuum")
    a("--pe", type=float, default=d.pe_bar, help="design exit pressure [bar] (default: = pamb)")
    a("--eps", type=float, default=d.eps, help="expansion ratio Ae/At (overrides --pe)")
    a("--ox", default=d.ox, help="oxidizer (RocketCEA name)")
    a("--fuel", default=d.fuel, help="fuel (RocketCEA name)")
    a("--of", type=float, default=d.of_ratio, help="core O/F mass ratio (film fuel not included)")
    a("--film", type=float, default=d.film_pct_core,
      help="fuel film flow as %% of the CORE mass flow (ox + core fuel), e.g. 15 = 0.15 kg/s per 1 kg/s core")
    a("--of-sweep", type=float, nargs=3, default=d.of_sweep, metavar=("MIN", "MAX", "STEP"),
      help="sweep O/F and tabulate")
    a("--optimize-of", action="store_true", default=d.optimize_of, help="design at the max-Isp O/F from the sweep")
    a("--frozen", action="store_true", default=d.frozen, help="frozen composition from the chamber")
    a("--frozen-at-throat", action="store_true", default=d.frozen_at_throat, help="equilibrium to throat, frozen after")
    a("--lstar", type=float, default=d.lstar_m, help="characteristic length L* [m]")
    a("--cr", type=float, default=d.contraction_ratio, help="contraction ratio Ac/At (default: Humble correlation)")
    a("--conv-angle", type=float, default=d.conv_half_angle_deg, help="convergent half angle [deg]")
    a("--cone-angle", type=float, default=d.cone_half_angle_deg, help="divergent cone half angle [deg]")
    a("--eta-cstar", type=float, default=d.eta_cstar, help="c* (injector/combustion) efficiency")
    a("--eta-nozzle", type=float, default=d.eta_nozzle, help="nozzle efficiency excl. divergence")
    a("--eta-div", type=float, default=d.eta_divergence, help="divergence efficiency override")
    a("--csv", default=d.csv_path, help="contour CSV output path ('' to skip)")
    a("--plot", action="store_true", default=d.plot, help="save a contour plot")
    a("--plot-path", default=d.plot_path)
    a("--show", action="store_true", help="show the plot interactively")
    o = ap.parse_args(argv)

    inp = EngineInputs(
        pc_bar=o.pc, thrust_n=o.thrust, pamb_bar=o.pamb, pe_bar=o.pe, eps=o.eps,
        ox=o.ox, fuel=o.fuel, of_ratio=o.of, film_pct_core=o.film,
        of_sweep=tuple(o.of_sweep) if o.of_sweep else None, optimize_of=o.optimize_of,
        frozen=o.frozen, frozen_at_throat=o.frozen_at_throat,
        lstar_m=o.lstar, contraction_ratio=o.cr, conv_half_angle_deg=o.conv_angle,
        cone_half_angle_deg=o.cone_angle,
        eta_cstar=o.eta_cstar, eta_nozzle=o.eta_nozzle, eta_divergence=o.eta_div,
        csv_path=o.csv or None, plot=o.plot or o.show, plot_path=o.plot_path,
    )
    return inp, o.show


def main(argv=None) -> None:
    inp, show = parse_args(argv)
    design = design_engine(inp)
    print_report(design)
    if inp.csv_path:
        write_csv(design, inp.csv_path)
    if inp.plot:
        plot_contour(design, inp.plot_path, show)


if __name__ == "__main__":
    main()
