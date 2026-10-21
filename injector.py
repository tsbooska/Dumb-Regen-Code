"""
Fuel-side pressure budget and plain-orifice injector sizing.

The fuel is the regen coolant: it flows through the cooling channels, into the injector
manifold, then through the fuel orifices into the chamber.

    Fuel:  P_regen_in --[regen channel friction]--> P_inj,f --[s_f * Pc]--> Pc
    Ox:                                             P_inj,ox --[s_ox * Pc]--> Pc

    P_inj = Pc (1 + s)            s = injector stiffness dP_inj / Pc (input per side)
    P_regen_in = P_inj,f + dP_regen   (solved by fixed-point around run_thermal, since
                                       coolant properties depend on pressure)

Orifices (plain holes, incompressible):
    A_tot = mdot / (Cd sqrt(2 rho dP_inj)),   d = sqrt(4 A_tot / (pi N)),   V = mdot / (rho A_tot)
    Fuel density at the regen outlet state (hot fuel); oxidizer density from CoolProp at the
    oxidizer inlet temperature and manifold pressure, or a user-supplied value.

Feed-line and pump losses upstream of the regen inlet are not included.

Usage
-----
    py -3.14 injector.py                       # engine + thermal + injector defaults
    py -3.14 injector.py --pc 30 --of 2.3      # any shitty_rocket_code.py flag is accepted
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

from CoolProp.CoolProp import PropsSI

import regen_thermal as rt
import wall_stress as ws
from shitty_rocket_code import BAR, EngineDesign, design_engine, parse_args, print_report


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------
@dataclass
class InjectorInputs:
    stiff_fuel: float = 0.20             # fuel injector dP / Pc
    stiff_ox: float = 0.20               # oxidizer injector dP / Pc
    n_fuel: int = 24                     # number of fuel orifices
    n_ox: int = 24                       # number of oxidizer orifices
    cd_fuel: float = 0.65                # discharge coefficient (sharp-edged ~0.6-0.65, chamfered ~0.8)
    cd_ox: float = 0.65
    t_ox_K: float = 90.0                 # oxidizer temperature at the injector manifold
    ox_density: float | None = None      # kg/m^3; None -> CoolProp from the oxidizer name
    n_film: int = 24                     # number of fuel film-cooling holes (used when the film is on)
    cd_film: float = 0.65


# CEA oxidizer name -> CoolProp fluid
OX_COOLPROP = {"LOX": "Oxygen", "GOX": "Oxygen", "O2": "Oxygen", "O2(L)": "Oxygen",
               "N2O": "NitrousOxide", "F2": "Fluorine", "F2(L)": "Fluorine"}


# --------------------------------------------------------------------------------------
# Pieces
# --------------------------------------------------------------------------------------
def ox_density(inj: InjectorInputs, ox: str, p: float) -> tuple[float, str]:
    if inj.ox_density is not None:
        return inj.ox_density, "user input"
    if ox not in OX_COOLPROP:
        raise ValueError(f"No CoolProp model for oxidizer '{ox}': enter its density at the injector.")
    name = OX_COOLPROP[ox]
    try:
        rho = PropsSI("D", "T", inj.t_ox_K, "P", p, name)
    except ValueError:  # exactly on the saturation line: take the liquid side
        rho = PropsSI("D", "T", inj.t_ox_K - 0.05, "P", p, name)
    return rho, f"CoolProp {name} at {inj.t_ox_K:.0f} K, {p / BAR:.1f} bar"


def size_orifices(mdot: float, rho: float, dp: float, cd: float, n: int) -> dict:
    area = mdot / (cd * math.sqrt(2 * rho * dp))
    return dict(area=area, d=math.sqrt(4 * area / (math.pi * n)), vel=mdot / (rho * area),
                mdot_hole=mdot / n, rho=rho, dp=dp, cd=cd, n=n, mdot=mdot)


def solve_feed_system(design: EngineDesign, th: rt.ThermalInputs, inj: InjectorInputs,
                      cooling: bool = True, stress: ws.StressInputs | None = None) -> dict:
    """Pressure budget for both sides, regen inlet pressure, and both orifice sets.

    With cooling=False the regen jacket is skipped (dP_regen = 0, fuel at its inlet temperature).
    With `stress` given (and cooling on), the hot-wall stresses of the converged solution are added.
    """
    warnings: list[str] = []
    pc = design.inputs.pc_bar * BAR
    p_inj_f = pc * (1 + inj.stiff_fuel)
    p_inj_o = pc * (1 + inj.stiff_ox)
    mdot_f, mdot_o = design.sizing["mdot_f"], design.sizing["mdot_ox"]

    res, iterations = None, 0
    film_on = design.sizing.get("mdot_film", 0.0) > 0
    if cooling:
        # Secant on the regen inlet pressure so the fuel arrives at the manifold at P_inj,f.
        # First step is the fixed-point update (exact if dP_regen didn't depend on pressure).
        # With film, the film leaves the fuel manifold at the regen outlet temperature, so each pass
        # also feeds the previous outlet temperature back as the film injection temperature.
        film_t = [th.t_in_K]

        def outlet_error(p_in):
            th.p_in_bar = p_in / BAR
            r = rt.run_thermal(design, th, film_t[0] if film_on else None)
            film_t.append(r["t_out"])
            return r["p_out"] - p_inj_f, r

        p_a = 2.0 * p_inj_f
        while True:  # make sure the first guess doesn't run out of pressure in the channels
            try:
                e_a, res = outlet_error(p_a)
                break
            except ValueError as exc:
                if "pressure dropped below zero" not in str(exc) or p_a > 100 * p_inj_f:
                    raise
                p_a *= 2.0
        p_b = p_a - e_a
        for iterations in range(2, 31):
            film_t[0] = film_t[-1]
            e_b, res = outlet_error(p_b)
            film_ok = not film_on or abs(film_t[-1] - film_t[0]) < 0.5   # film inlet T converged (K)
            if abs(e_b) < 100.0 and film_ok:  # 1e-3 bar
                break
            slope = (e_b - e_a) / (p_b - p_a) if p_b != p_a else 1.0
            p_a, e_a = p_b, e_b
            p_b -= e_b / (slope if abs(slope) > 1e-3 else 1.0)
        else:
            warnings.append("Regen inlet pressure did not converge to 1e-3 bar.")
        cool, t_f, dp_regen = res["cool"], res["t_out"], res["dp_regen"]
    else:
        cool, t_f, dp_regen = rt.make_coolant(design, th, warnings), th.t_in_K, 0.0

    rho_f = cool.props(t_f, p_inj_f)["rho"]
    rho_o, ox_src = ox_density(inj, design.inputs.ox, p_inj_o)
    fuel = size_orifices(mdot_f, rho_f, p_inj_f - pc, inj.cd_fuel, inj.n_fuel)
    fuel.update(t=t_f, source=f"{cool.name} at {t_f:.0f} K, {p_inj_f / BAR:.1f} bar"
                + (" (regen outlet)" if cooling else " (no regen)"))
    ox = size_orifices(mdot_o, rho_o, p_inj_o - pc, inj.cd_ox, inj.n_ox)
    ox.update(t=inj.t_ox_K, source=ox_src)
    film = None
    if film_on:
        # Film holes are fed from the fuel manifold: same pressure drop and fuel density as the fuel orifices
        film = size_orifices(design.sizing["mdot_film"], rho_f, p_inj_f - pc, inj.cd_film, inj.n_film)
        film.update(t=t_f, source=fuel["source"])

    for name, s in (("Fuel", inj.stiff_fuel), ("Oxidizer", inj.stiff_ox)):
        if s < 0.15:
            warnings.append(f"{name} stiffness {s:.0%} < 15%: weak isolation from chamber pressure "
                            "oscillations (chugging risk).")
        elif s > 0.40:
            warnings.append(f"{name} stiffness {s:.0%} > 40%: costs a lot of feed pressure for little gain.")
    for name, o in (("Fuel", fuel), ("Oxidizer", ox)) + ((("Film", film),) if film else ()):
        if o["d"] < 0.3e-3:
            warnings.append(f"{name} orifices are {o['d'] * 1e3:.2f} mm: hard to drill and easy to clog. "
                            "Use fewer holes.")
    if rho_o < 100:
        warnings.append(f"Oxidizer density {rho_o:.1f} kg/m^3 means a gas: the incompressible orifice "
                        "model does not apply (needs compressible / choked-flow sizing).")
    t_sat = cool.t_sat(p_inj_f)
    if t_sat is not None and t_f > t_sat - 10:
        warnings.append(f"Fuel reaches the injector at {t_f:.0f} K, within 10 K of saturation "
                        f"({t_sat:.0f} K at {p_inj_f / BAR:.1f} bar): flashing in the orifices likely.")

    budget = dict(pc=pc, p_inj_f=p_inj_f, p_inj_o=p_inj_o, dp_inj_f=p_inj_f - pc, dp_inj_o=p_inj_o - pc,
                  dp_regen=dp_regen, p_feed_f=p_inj_f + dp_regen, p_feed_o=p_inj_o, iterations=iterations)
    momentum_ratio = (mdot_o * ox["vel"]) / (mdot_f * fuel["vel"])
    stresses = ws.wall_stresses(design, th, res, stress) if cooling and stress is not None else None
    return dict(thermal=res, budget=budget, fuel=fuel, ox=ox, film=film, momentum_ratio=momentum_ratio,
                warnings=warnings, inputs=inj, cooling=cooling, stress=stresses)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def print_feed_report(design: EngineDesign, feed: dict) -> None:
    b, f, o, inj = feed["budget"], feed["fuel"], feed["ox"], feed["inputs"]
    i = design.inputs
    bar = 1 / BAR
    print("=" * 70)
    print(f" FEED SYSTEM & INJECTOR: {i.fuel} (through regen) / {i.ox}")
    print("=" * 70)
    print(f"  {'PRESSURE BUDGET':<28}{'fuel':>12}{'ox':>12}")
    print(f"  {'Chamber pressure Pc':<28}{b['pc'] * bar:12.2f}{b['pc'] * bar:12.2f}  bar")
    print(f"  {'Injector dP':<28}{b['dp_inj_f'] * bar:12.2f}{b['dp_inj_o'] * bar:12.2f}  bar"
          f"   ({inj.stiff_fuel:.0%} / {inj.stiff_ox:.0%} of Pc)")
    print(f"  {'Injector manifold P':<28}{b['p_inj_f'] * bar:12.2f}{b['p_inj_o'] * bar:12.2f}  bar")
    regen = f"{b['dp_regen'] * bar:12.2f}" if feed["cooling"] else f"{'off':>12}"
    print(f"  {'Regen channel dP':<28}{regen}{'-':>12}  bar")
    print(f"  {'REQUIRED FEED PRESSURE':<28}{b['p_feed_f'] * bar:12.2f}{b['p_feed_o'] * bar:12.2f}  bar"
          "   (fuel at regen inlet, ox at manifold)")
    if feed["cooling"]:
        print(f"  {'':<28}(regen inlet pressure converged in {b['iterations']} passes)")

    sets = [f, o] + ([feed["film"]] if feed.get("film") else [])

    def row(label, fmt, key, scale=None, unit=""):
        vals = (s[key] if scale is None else s[key] * scale for s in sets)
        print(f"  {label:<28}" + "".join(format(v, fmt) for v in vals) + (f"  {unit}" if unit else ""))

    print(f"\n  {'ORIFICES (plain holes)':<28}{'fuel':>12}{'ox':>12}" + (f"{'film':>12}" if len(sets) > 2 else ""))
    row("Number of holes", "12d", "n")
    row("Discharge coefficient Cd", "12.3f", "cd")
    row("Mass flow", "12.4f", "mdot", unit="kg/s")
    row("Flow per hole", "12.3f", "mdot_hole", 1e3, "g/s")
    row("Density at injector", "12.1f", "rho", unit="kg/m^3")
    row("HOLE DIAMETER", "12.3f", "d", 1e3, "mm")
    row("Total orifice area", "12.2f", "area", 1e6, "mm^2")
    row("Velocity in hole", "12.2f", "vel", unit="m/s")
    print(f"  {'Momentum ratio ox/fuel':<28}{feed['momentum_ratio']:12.2f}" + ("   (core flows)" if len(sets) > 2 else ""))
    print(f"\n  Fuel density: {f['source']}")
    print(f"  Ox density:   {o['source']}")
    if len(sets) > 2:
        print(f"  Film holes:   fed from the fuel manifold (same dP and density as the fuel orifices); "
              f"film flow = {design.inputs.film_pct_core:g} % of core flow")
    if feed["warnings"]:
        print("\n WARNINGS")
        for msg in feed["warnings"]:
            print(f"  ! {msg}")
    print()


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)
    ap.add_argument("--channels", help="channel table CSV: x_mm, width_mm, height_mm, wall_t_mm")
    ap.add_argument("--coolant-table", help="fuel property CSV: T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK")
    o, engine_argv = ap.parse_known_args(argv)

    inp, _ = parse_args(engine_argv)
    th = rt.ThermalInputs(csv_path=None, plot=False, coolant_table=o.coolant_table)
    if o.channels:
        th.channel_table = rt.read_channel_csv(o.channels)

    design = design_engine(inp)
    print_report(design)
    feed = solve_feed_system(design, th, InjectorInputs())
    rt.print_thermal_report(design, th, feed["thermal"])
    print_feed_report(design, feed)


if __name__ == "__main__":
    main()
