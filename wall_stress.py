"""
Hot-wall (liner) stresses, safety factors, low-cycle fatigue life and buckling margin for the regen
jacket, station by station.

Loads, after Huzel & Huang, NASA SP-125 (2nd ed., 1971), Eqs. (4-27), (4-28), (4-31):

    sigma_b = (dp / 2) (w / t_w)^2              pressure bending, wall = beam built in at both ribs
                                                (Roark: root moment dp w^2/12, mid-span dp w^2/24,
                                                sigma = 6M/t^2; mid-span stress is half the root value)
    sigma_t = E alpha q t_w / (2 (1 - nu) k)    through-wall thermal stress (linear gradient, plate
            = E alpha (T_wg - T_wc) / (2 (1 - nu))   restrained against bending); q t_w / k = T_wg - T_wc
    sigma_l = E alpha (T_wg - T_jacket)         longitudinal thermal stress, gas-side wall restrained by
                                                the jacket at its outer temperature T_jacket (default
                                                ambient); compressive in the hot wall. alpha is the mean
                                                CTE from room temperature, taken at T_wg

    dp = p_coolant(x) - p_gas(x),   p_gas = Pc (1 + (g-1)/2 M^2)^(-g/(g-1))  (local static pressure)
    E, alpha, nu at the mean wall temperature; strengths at the temperature of the face checked.

Signed stresses are evaluated at four points per station: hot face and cold face, each at the rib root
and at channel mid-span. With dp > 0 pushing the wall toward the gas, bending puts the cold face in
tension at the rib root and the hot face in tension at mid-span. sigma_t is compressive on the hot face
and tensile on the cold face; sigma_l is compressive. At each point, with sigma_z = 0 at the free face
and no shear, the principal stresses are sigma_x, sigma_y and 0:

    sigma_x  = tangential stress (bending + through-wall thermal, signed for the point)
    sigma_y  = longitudinal stress = -sigma_l
    sigma_vm = sqrt(sigma_x^2 - sigma_x sigma_y + sigma_y^2)        (von Mises, = thesis Eq. 10 with signs)
    sigma_1  = max(sigma_x, sigma_y, 0) (most tensile),  sigma_3 = min(sigma_x, sigma_y, 0) (most compressive)
    sigma_p  = maximum principal stress = the principal stress of largest magnitude (sigma_1 or sigma_3)
    SF_yield = sigma_yield(T_face) / sigma_vm
No ultimate-strength factor is reported: the elastic stresses exceed yield in a regen liner, where they
are not real stresses, so a ratio to the ultimate strength would not mean anything.

The thesis magnitude sum, sqrt(0.5[(|sb|+st-sl)^2 + sl^2 + (|sb|+st)^2]), is reported for reference.

Low-cycle fatigue: each firing is one cycle from rest to steady firing and back. The total strain range
at each point is the elastic equivalent strain d_eps = sigma_vm / E, and the life N comes from the
material's measured strain-life curve, or from Manson's universal slopes where only tensile data exist
(materials.cycles_to_failure). Allowable firings = N / life_factor.

Longitudinal buckling of the gas-side wall, SP-125 Eq. (4-29), printed there with a typo in the
denominator; the reduced-modulus form (Gerard) is
    S_c = 4 E_t E_c t / ((sqrt(E_t) + sqrt(E_c))^2 sqrt(3 (1 - nu^2)) r)
Without tangent moduli, E_t = E_c = E, which gives classical elastic axial buckling of a cylinder,
S_c = E t / (r sqrt(3 (1 - nu^2))), with r the mean wall radius. SP-125 keeps sigma_l <= 0.9 S_c, so the
margin is 0.9 S_c / sigma_l - 1. Where sigma_l exceeds yield the tangent modulus is far below E, so the
elastic margin is an upper bound (unconservative); this is flagged.

Elastic analysis only: copper liners normally yield in compression on the hot face each firing; the
fatigue estimate uses elastic strains and no ratcheting or creep.

Usage
-----
    py -3.14 wall_stress.py                          # engine, thermal, injector and stress defaults
    py -3.14 wall_stress.py --material GRCop-42      # any shitty_rocket_code.py flag too
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

import materials as mat
from shitty_rocket_code import BAR, EngineDesign

# (face, position, sign of the bending stress per unit root stress, sign of the through-wall thermal stress)
POINTS = (("hot face", "rib root", -1.0, -1.0), ("hot face", "mid-span", 0.5, -1.0),
          ("cold face", "rib root", 1.0, 1.0), ("cold face", "mid-span", -0.5, 1.0))
COMPONENTS = ("tangential", "longitudinal")

ELASTIC_DISCLAIMER = (
    "ALL STRUCTURAL RESULTS ASSUME FULLY ELASTIC BEHAVIOUR (sigma = E epsilon), even past yield. Stresses "
    "above the yield strength are elastic-equivalent values, not real stresses: the actual stress stays near "
    "yield and the excess becomes plastic strain. A yield safety factor below 1 means the wall yields, not "
    "how much. Fatigue life uses the elastic strain range (close to the true total strain for these "
    "thermally driven loads, but without plastic strain concentration, ratcheting or creep). The buckling "
    "margin uses E in place of the tangent modulus and is optimistic wherever the wall has yielded.")

# Explanations printed with the report (and shown in the GUI's Structures tab)
DEFINITIONS = [
    ("dp", "Pressure difference across the hot wall: coolant pressure minus the local static pressure of the "
           "gas, p_gas = Pc (1 + (g-1)/2 M^2)^(-g/(g-1)). Positive pushes the wall toward the gas."),
    ("s_b,root / s_b,mid", "Pressure bending stress. The wall between two ribs is a beam built in at both "
           "ends, span w (channel width), depth t_w, uniform load dp. Rib root: s_b = (dp/2)(w/t_w)^2 (thesis "
           "Eq. 6), tensile on the cold face. Mid-span: half that, tensile on the hot face. (Roark, Table 8.1.)"),
    ("s_t", "Through-wall thermal stress: s_t = E a (T_wg - T_wc) / (2(1 - nu)) = E a q t_w / (2(1 - nu) k) "
           "(thesis Eq. 7; SP-125 Eqs. 4-27, 4-31). Linear temperature drop across a plate held flat. "
           "Compressive on the hot face, tensile on the cold face, acting in the tangential direction."),
    ("s_l", "Longitudinal thermal stress: s_l = E a (T_wg - T_jacket) (thesis Eq. 8; SP-125 Eq. 4-28). The "
           "gas-side wall is held axially by the jacket at its outer temperature T_jacket. Compressive."),
    ("s_x", "Tangential stress at a check point = signed s_b at that point + signed s_t."),
    ("s_y", "Longitudinal stress at a check point = -s_l (the same at all four points)."),
    ("s_vm", "von Mises equivalent stress, sqrt(s_x^2 - s_x s_y + s_y^2) with s_z = 0 at the free face "
           "(thesis Eqs. 9-10, with signs). Compared with the yield strength."),
    ("s_1 / s_3", "Principal stresses: s_1 = most tensile = max(s_x, s_y, 0); s_3 = most compressive = "
           "min(s_x, s_y, 0). With no shear on the faces, s_x and s_y are principal stresses."),
    ("s_p", "Maximum principal stress = the principal stress of largest magnitude (s_1 or s_3, keeping its "
           "sign), with its direction (tangential or longitudinal)."),
    ("SF_y", "Yield safety factor = yield strength at the face temperature / s_vm."),
    ("strain", "Total strain range per firing, elastic equivalent: s_vm / E (rest -> firing -> rest)."),
    ("N_f", "Low-cycle fatigue life (cycles to failure) from the material's strain-life curve, or Manson's "
           "universal slopes (estimate). Allowable firings = N_f / life factor."),
    ("S_c / MS_b", "Longitudinal buckling of the gas-side wall (SP-125 Eq. 4-29 with E_t = E_c = E): "
           "S_c = E t_w / (r sqrt(3(1 - nu^2))); margin MS_b = 0.9 S_c / s_l - 1 (SP-125 limit s_l <= 0.9 S_c). "
           "Elastic, so optimistic wherever s_l exceeds yield."),
]


@dataclass
class StressInputs:
    material: str = "CuCrZr"            # materials.py name or GUI label
    sf_yield_target: float = 1.0         # for warnings and summary colours
    t_jacket_K: float = 293.15           # jacket outer temperature for sigma_l; default ambient
    life_factor: float = 4.0             # divides the fatigue life to give allowable firings


def _region(design: EngineDesign, x: float) -> str:
    if abs(x) < 0.02 * design.sizing["rt"]:
        return "throat"
    if x > 0:
        return "nozzle"
    return "convergent" if x > -design.chamber["l_conv"] else "chamber"


def wall_stresses(design: EngineDesign, th, res: dict, si: StressInputs) -> dict:
    """Stresses, safety factors, fatigue life and buckling margin at every cooled station.

    Per station (arrays of length n): dp, p_gas, sb (root, signed), sb_mid, st, sl, E, alpha, nu, t_mean,
    sc, bm, thesis_vm, and the station envelope vm (max), sp (largest |s_p|), sfy, nf (min), eps.
    Per point (arrays n x 4, points in POINTS order): sx, sy, vm_p, s1_p, s3_p, sp_p, sfy_p, eps_p,
    nf_p, t_face.
    """
    warnings: list[str] = []
    m = mat.get(si.material)
    st, o, cooled = res["st"], res["out"], res["cooled"]
    n = len(st["x"])
    g = design.perf["gam_t"]
    pc = design.inputs.pc_bar * BAR
    p_gas = pc * (1 + 0.5 * (g - 1) * st["mach"] ** 2) ** (-g / (g - 1))

    if abs(th.k_wall - m.k) > 0.05 * m.k:
        warnings.append(f"Thermal k_wall = {th.k_wall:g} W/m-K differs from the {m.name} room-temperature value "
                        f"({m.k:.0f}). Stresses use the temperatures from the thermal solution as computed.")

    out = {k: np.full(n, np.nan) for k in ("dp", "sb", "sb_mid", "st", "sl", "E", "alpha", "nu", "t_mean",
                                           "sc", "bm", "thesis_vm", "vm", "sp", "sfy", "nf", "eps")}
    pts = {k: np.full((n, 4), np.nan) for k in ("sx", "sy", "vm_p", "s1_p", "s3_p", "sp_p", "sfy_p",
                                                "eps_p", "nf_p", "t_face")}
    sp_comp = np.zeros((n, 4), dtype=int)       # 0 tangential, 1 longitudinal
    life_method = "no fatigue data"
    inelastic = False

    for j in np.flatnonzero(cooled):
        w, t = res["w"][j], res["t"][j]
        t_wg, t_wc = o["t_wg"][j], o["t_wc"][j]
        t_m = 0.5 * (t_wg + t_wc)
        e, a, nu = (mat.prop(m.name, k, t_m, warnings) for k in ("E", "alpha", "nu"))
        dp = o["p_b"][j] - p_gas[j]
        sb = 0.5 * dp * (w / t) ** 2
        s_t = e * a * (t_wg - t_wc) / (2 * (1 - nu))
        s_l = e * mat.prop(m.name, "alpha", t_wg, warnings) * (t_wg - si.t_jacket_K)
        for k, v in (("dp", dp), ("sb", sb), ("sb_mid", 0.5 * sb), ("st", s_t), ("sl", s_l), ("E", e),
                     ("alpha", a), ("nu", nu), ("t_mean", t_m)):
            out[k][j] = v
        b = abs(sb)
        out["thesis_vm"][j] = np.sqrt(0.5 * ((b + s_t - s_l) ** 2 + s_l ** 2 + (b + s_t) ** 2))

        # Longitudinal buckling of the gas-side wall (SP-125 Eq. 4-29 with E_t = E_c = E)
        r_m = st["r"][j] + 0.5 * t
        sc = e * t / (r_m * np.sqrt(3 * (1 - nu * nu)))
        out["sc"][j], out["bm"][j] = sc, 0.9 * sc / s_l - 1 if s_l > 0 else np.inf
        ys_hot = mat.prop(m.name, "sy", t_wg, warnings)
        inelastic |= ys_hot is not None and s_l > ys_hot

        for p, (face, _, kb, kt) in enumerate(POINTS):
            sx, sy = kb * sb + kt * s_t, -s_l
            vm = np.sqrt(sx * sx - sx * sy + sy * sy)
            s1, s3 = max(sx, sy, 0.0), min(sx, sy, 0.0)
            sp_comp[j, p] = 0 if abs(sx) >= abs(sy) else 1
            sp = sx if sp_comp[j, p] == 0 else sy
            t_face = t_wg if face == "hot face" else t_wc
            ys = mat.prop(m.name, "sy", t_face, warnings)
            life, life_method = mat.cycles_to_failure(m.name, vm / e, t_face, warnings)
            for k, v in (("sx", sx), ("sy", sy), ("vm_p", vm), ("s1_p", s1), ("s3_p", s3), ("sp_p", sp),
                         ("sfy_p", ys / vm if ys is not None and vm > 0 else np.nan), ("eps_p", vm / e), ("nf_p", np.nan if life is None else life), ("t_face", t_face)):
                pts[k][j, p] = v

    # Station envelopes and the point that sets them
    env = {}
    with np.errstate(all="ignore"):
        for key, arr, pick in (("vm", pts["vm_p"], np.nanargmax), ("sp", np.abs(pts["sp_p"]), np.nanargmax),
                               ("sfy", pts["sfy_p"], np.nanargmin), ("nf", pts["nf_p"], np.nanargmin)):
            idx = np.full(n, -1, dtype=int)
            for j in range(n):
                if np.isfinite(arr[j]).any():
                    idx[j] = int(pick(arr[j]))
            env[key] = idx
    for key, src in (("vm", "vm_p"), ("sp", "sp_p"), ("sfy", "sfy_p"), ("nf", "nf_p")):
        ok = env[key] >= 0
        out[key][ok] = pts[src][ok, env[key][ok]]
    ok = env["nf"] >= 0
    out["eps"][ok] = pts["eps_p"][ok, env["nf"][ok]]

    def where(j: int, p: int, with_component: bool = False) -> dict:
        face, pos = POINTS[p][:2]
        loc = dict(j=j, p=p, x=float(st["x"][j]), region=_region(design, st["x"][j]), face=face, position=pos,
                   t_face=float(pts["t_face"][j, p]))
        text = f"{face} at {pos}, {loc['region']}, x = {loc['x'] * 1e3:.1f} mm"
        if with_component:
            sp = pts["sp_p"][j, p]
            loc["component"] = COMPONENTS[sp_comp[j, p]]
            loc["sense"] = "tensile" if sp > 0 else "compressive"
            text = f"{loc['component']}, {loc['sense']}, " + text
        loc["text"] = text
        return loc

    extremes = {}
    for key, pick, comp in (("vm", np.nanargmax, False), ("sfy", np.nanargmin, False), ("nf", np.nanargmin, False)):
        if np.isfinite(out[key]).any():
            j = int(pick(out[key]))
            extremes[key] = dict(value=out[key][j], **where(j, env[key][j], comp))
    if np.isfinite(out["sp"]).any():
        j = int(np.nanargmax(np.abs(out["sp"])))
        extremes["sp"] = dict(value=out["sp"][j], **where(j, env["sp"][j], True))
    if "nf" in extremes:
        extremes["nf"]["eps"] = out["eps"][extremes["nf"]["j"]]
    if np.isfinite(out["thesis_vm"]).any():
        j = int(np.nanargmax(out["thesis_vm"]))
        extremes["thesis_vm"] = dict(value=out["thesis_vm"][j], j=j, x=float(st["x"][j]),
                                     region=_region(design, st["x"][j]))
    if np.isfinite(out["bm"]).any():
        j = int(np.nanargmin(out["bm"]))
        extremes["bm"] = dict(value=out["bm"][j], sc=out["sc"][j], sl=out["sl"][j], j=j, x=float(st["x"][j]),
                              region=_region(design, st["x"][j]),
                              text=f"gas-side wall, {_region(design, st['x'][j])}, x = {st['x'][j] * 1e3:.1f} mm")

    if "nf" in extremes and extremes["nf"]["value"] / si.life_factor < 1:
        warnings.append(f"Fatigue life {extremes['nf']['value']:.0f} cycles is below the life factor "
                        f"({si.life_factor:g}): fewer than one allowable firing ({extremes['nf']['text']}).")
    if life_method == "no fatigue data":
        warnings.append(f"{m.name}: no fatigue data (no strain-life curve or reduction of area); life not computed.")
    if "bm" in extremes and extremes["bm"]["value"] < 0:
        warnings.append(f"Longitudinal buckling margin {extremes['bm']['value']:.2f} < 0 ({extremes['bm']['text']}).")
    if inelastic:
        warnings.append("Longitudinal stress exceeds yield on the hot wall: inelastic buckling governs and the "
                        "elastic buckling margin (E_t = E_c = E) is unconservative there.")
    if "sfy" in extremes and extremes["sfy"]["value"] < si.sf_yield_target:
        warnings.append(f"Yield safety factor {extremes['sfy']['value']:.2f} < {si.sf_yield_target:g} "
                        f"({extremes['sfy']['text']}). Elastic analysis: the wall yields locally; see the "
                        "fatigue life.")
    if m.sy is None:
        warnings.append(f"{m.name}: no yield-strength data, yield safety factor not computed.")

    return dict(material=m, inputs=si, p_gas=p_gas, envelope_point=env, sp_component=sp_comp, extremes=extremes,
                warnings=warnings, life_method=life_method, **out, **pts)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def _wrap(text: str, width: int = 96, indent: str = "      ") -> str:
    words, lines, line = text.split(), [], ""
    for wd in words:
        if len(line) + len(wd) + 1 > width:
            lines.append(line)
            line = wd
        else:
            line = f"{line} {wd}".strip()
    lines.append(line)
    return ("\n" + indent).join(lines)


def print_stress_report(design: EngineDesign, res: dict, sres: dict, definitions: bool = True) -> None:
    m, ex, si, x = sres["material"], sres["extremes"], sres["inputs"], res["st"]["x"]
    o, mp = res["out"], 1e-6
    print("=" * 100)
    print(f" STRUCTURAL ANALYSIS OF THE HOT WALL: {m.name}")
    print("=" * 100)
    print(f"  {_wrap(ELASTIC_DISCLAIMER, 94, '  ')}")
    print("\n 1. INPUTS AND LOAD CASE")
    print(f"  Material                     {m.name}, {m.condition}")
    print(f"  Load case                    steady firing; dp = coolant - local gas static pressure")
    print(f"  Jacket outer temperature     {si.t_jacket_K:.1f} K   (restrains the gas-side wall axially, s_l)")
    print(f"  Yield safety-factor target   {si.sf_yield_target:g}   (vs von Mises)")
    print(f"  Fatigue life factor          {si.life_factor:g}   (allowable firings = N_f / factor)")
    print(f"  Fatigue method               {sres['life_method']}")
    print("  Check points                 hot face and cold face, each at the rib root and at mid-span")
    print("  Sign convention              positive = tension; x = tangential, y = longitudinal, z = radial (= 0)")

    if definitions:
        print("\n 2. DEFINITIONS OF THE CHECKED STRESSES AND OUTPUTS")
        for sym, text in DEFINITIONS:
            print(f"  {sym:<20}{_wrap(text, 78, ' ' * 22)}")

    rows = [("injector", 0), ("conv. start", int(np.argmin(abs(x + design.chamber["l_conv"])))),
            ("throat", int(np.argmin(abs(x)))), ("exit", len(x) - 1)]
    for name, key in (("max s_vm", "vm"), ("max |s_p|", "sp"), ("min SF_y", "sfy"), ("min N_f", "nf"),
                      ("min MS_b", "bm")):
        if key in ex and ex[key]["j"] not in [j for _, j in rows]:
            rows.append((name, ex[key]["j"]))

    print("\n 3. LOAD COMPONENTS AT KEY STATIONS (per station, before combining)")
    print(f"  {'station':<12}{'x':>7}{'T_wg':>7}{'T_wc':>7}{'dp':>7}{'s_b,root':>10}{'s_b,mid':>9}{'s_t':>8}"
          f"{'s_l':>8}{'E':>7}{'alpha':>8}{'nu':>6}{'S_c':>8}{'MS_b':>7}")
    print(f"  {'':<12}{'mm':>7}{'K':>7}{'K':>7}{'bar':>7}{'MPa':>10}{'MPa':>9}{'MPa':>8}{'MPa':>8}{'GPa':>7}"
          f"{'1e-6/K':>8}{'':>6}{'MPa':>8}{'':>7}")
    for name, j in rows:
        print(f"  {name:<12}{x[j] * 1e3:7.1f}{o['t_wg'][j]:7.0f}{o['t_wc'][j]:7.0f}{sres['dp'][j] / BAR:7.2f}"
              f"{sres['sb'][j] * mp:10.2f}{sres['sb_mid'][j] * mp:9.2f}{sres['st'][j] * mp:8.1f}"
              f"{sres['sl'][j] * mp:8.1f}{sres['E'][j] / 1e9:7.1f}{sres['alpha'][j] * 1e6:8.2f}{sres['nu'][j]:6.3f}"
              f"{sres['sc'][j] * mp:8.0f}{sres['bm'][j]:7.2f}")

    print("\n 4. STRESSES AT EACH CHECK POINT, KEY STATIONS")
    print(f"  {'station':<12}{'point':<21}{'T':>6}{'s_x':>8}{'s_y':>8}{'s_vm':>8}{'s_1':>8}{'s_3':>8}{'s_p':>8}"
          f"{'SF_y':>7}{'strain':>8}{'N_f':>9}")
    print(f"  {'':<12}{'':<21}{'K':>6}{'MPa':>8}{'MPa':>8}{'MPa':>8}{'MPa':>8}{'MPa':>8}{'MPa':>8}{'':>7}"
          f"{'%':>8}{'cycles':>9}")
    for name, j in rows:
        for p, (face, pos, _, _) in enumerate(POINTS):
            label = name if p == 0 else ""
            print(f"  {label:<12}{face + ', ' + pos:<21}{sres['t_face'][j, p]:6.0f}{sres['sx'][j, p] * mp:8.1f}"
                  f"{sres['sy'][j, p] * mp:8.1f}{sres['vm_p'][j, p] * mp:8.1f}{sres['s1_p'][j, p] * mp:8.1f}"
                  f"{sres['s3_p'][j, p] * mp:8.1f}{sres['sp_p'][j, p] * mp:8.1f}{sres['sfy_p'][j, p]:7.2f}"
                  f"{sres['eps_p'][j, p] * 100:8.3f}{sres['nf_p'][j, p]:9.3g}")

    print("\n 5. GOVERNING VALUES (whole jacket)")
    for key, label, unit in (("vm", "Max von Mises stress s_vm", "MPa"),
                             ("sp", "Max principal stress s_p", "MPa"),
                             ("sfy", "Min yield SF (vs s_vm)", "")):
        if key in ex:
            v = ex[key]["value"] * (mp if unit else 1)
            print(f"  {label:<30}{v:9.2f} {unit:<4} {ex[key]['text']}, T = {ex[key]['t_face']:.0f} K")
    if "thesis_vm" in ex:
        print(f"  {'Thesis Eq. 10 (magnitudes)':<30}{ex['thesis_vm']['value'] * mp:9.2f} MPa  "
              f"{ex['thesis_vm']['region']}, x = {ex['thesis_vm']['x'] * 1e3:.1f} mm (reference only)")

    print("\n 6. LOW-CYCLE FATIGUE")
    if "nf" in ex:
        e = ex["nf"]
        print(f"  {'Min cycles to failure N_f':<30}{e['value']:9.3g}       {e['text']}")
        print(f"  {'Strain range there':<30}{e['eps'] * 100:9.3f} %     s_vm / E")
        print(f"  {'Allowable firings':<30}{e['value'] / si.life_factor:9.3g}       N_f / {si.life_factor:g}")
        print(f"  {'Method':<30}{sres['life_method']}: {m.lcf.source if m.lcf else mat.MANSON}")
    else:
        print(f"  {sres['life_method']}")

    print("\n 7. LONGITUDINAL BUCKLING")
    if "bm" in ex:
        e = ex["bm"]
        print(f"  {'Min margin MS_b':<30}{e['value']:9.2f}       {e['text']}")
        print(f"  {'s_l vs 0.9 S_c there':<30}{e['sl'] * mp:9.0f} MPa   vs {0.9 * e['sc'] * mp:.0f} MPa "
              "(elastic critical stress, E_t = E_c = E)")

    if sres["warnings"]:
        print("\n WARNINGS")
        for msg in sres["warnings"]:
            print(f"  ! {msg}")
    print()


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv=None) -> None:
    import injector as inj
    import regen_thermal as rt
    from shitty_rocket_code import design_engine, parse_args, print_report

    d = StressInputs()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)
    ap.add_argument("--material", default=d.material, help=f"wall material: {', '.join(mat.MATERIALS)}")
    ap.add_argument("--sf-yield", type=float, default=d.sf_yield_target, help="yield safety factor target")
    ap.add_argument("--t-jacket", type=float, default=d.t_jacket_K,
                    help="jacket outer temperature for the longitudinal stress [K]")
    ap.add_argument("--life-factor", type=float, default=d.life_factor,
                    help="fatigue life factor: allowable firings = cycles to failure / factor")
    ap.add_argument("--channels", help="channel table CSV: x_mm, width_mm, height_mm, wall_t_mm")
    ap.add_argument("--coolant-table", help="fuel property CSV: T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK")
    o, engine_argv = ap.parse_known_args(argv)

    inp, _ = parse_args(engine_argv)
    inp.csv_path = None
    th = rt.ThermalInputs(csv_path=None, plot=False, coolant_table=o.coolant_table)
    if o.channels:
        th.channel_table = rt.read_channel_csv(o.channels)
    si = StressInputs(material=o.material, sf_yield_target=o.sf_yield, t_jacket_K=o.t_jacket,
                      life_factor=o.life_factor)

    design = design_engine(inp)
    print_report(design)
    feed = inj.solve_feed_system(design, th, inj.InjectorInputs(), stress=si)
    rt.print_thermal_report(design, th, feed["thermal"])
    inj.print_feed_report(design, feed)
    print_stress_report(design, feed["thermal"], feed["stress"])


if __name__ == "__main__":
    main()
