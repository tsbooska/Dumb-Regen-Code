"""
Liner (hot-wall) material library for the regen jacket.

Each material carries its room-temperature thermal conductivity (the thermal circuit uses a single,
constant k) and, versus temperature, Young's modulus E, the mean (secant) coefficient of thermal
expansion alpha from room temperature, Poisson's ratio nu, the 0.2 % yield strength and the ultimate
tensile strength. Every property records its source. Values read from a figure, computed from a
fitted equation, or held constant for lack of temperature data say so in their source string.

Lookups interpolate linearly in temperature. Outside a table's range the end value is used and a
warning is added; a property tabulated at a single temperature is held constant (with a warning).

    prop("GRCop-42", "sy", 700.0, warnings)   # yield strength [Pa] at 700 K
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

C0 = 273.15           # K at 0 °C
MPA, GPA = 1e6, 1e9
KSI = 6.894757e6      # Pa
BTU_IN = 0.1442279    # W/m-K per Btu·in/(ft²·h·°F)


def _f2k(f: float) -> float:
    return (f - 32.0) / 1.8 + C0


@dataclass(frozen=True)
class Table:
    t: tuple[float, ...]     # K
    v: tuple[float, ...]     # SI
    source: str


def _tab(source: str, pairs, scale: float = 1.0, t_unit: str = "C") -> Table:
    conv = {"C": lambda x: x + C0, "F": _f2k, "K": lambda x: x}[t_unit]
    pairs = sorted(pairs)
    return Table(tuple(conv(t) for t, _ in pairs), tuple(v * scale for _, v in pairs), source)


def _secant(source: str, strain_pairs, t_ref_c: float) -> Table:
    """Mean CTE from t_ref from (T °C, thermal strain) pairs."""
    return _tab(source, [(t, e / (t - t_ref_c)) for t, e in strain_pairs if t > t_ref_c])


@dataclass(frozen=True)
class LCF:
    """Measured strain-life curve: total strain range = a * N^b (fully reversed, strain controlled)."""
    a: float
    b: float
    n_range: tuple[float, float]    # cycles covered by the test data
    t_max: float                    # K, highest test temperature
    source: str


@dataclass(frozen=True)
class Material:
    name: str
    label: str               # GUI label
    condition: str
    k: float                 # W/m-K at room temperature
    k_source: str
    E: Table                 # Pa
    alpha: Table             # 1/K, mean from room temperature
    nu: Table
    sy: Table | None         # Pa, 0.2 % yield
    su: Table | None         # Pa, ultimate tensile
    lcf: LCF | None = None   # measured strain-life curve
    ra: Table | None = None  # reduction of area [-], for the universal-slopes estimate when lcf is None


# --------------------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------------------
UKAEA = ("Zhang, Gaganidze, Gorley, 'Development of the material property handbook and database of "
         "CuCrZr', Fusion Eng. Des. 144 (2019) 148-153")
BAKER = ("Baker, Fiedler, Rosler, 'Stress evolution in thermal barrier coatings for rocket engine "
         "applications', Mech. Adv. Mater. Mod. Process. 1:5 (2015), Table 1 (CuCr1Zr; Kupferinstitut 2005, "
         "Dies 1967, Fassin 2015)")
GRADL = ("Gradl, Protz, Cooper, Garcia, Ellis, Evans, 'GRCop-42 Development and Hot-fire Testing Using "
         "Additive Manufacturing Powder Bed Fusion for Channel-Cooled Combustion Chambers', AIAA 2019-4228 "
         "(NTRS 20190030433)")
MAHTABI = ("Mahtabi et al., 'Multi-Scale Mechanical Characterization of Additively Manufactured GRCop-42 and "
           "GRCop-84 Alloys', JOM 77 (2025) 9973-9989")
NIST177 = ("Simon, Drexler, Reed, 'Properties of Copper and Copper Alloys at Cryogenic Temperatures', "
           "NIST Monograph 177 (1992)")
NASA_TB = "NASA Tech Brief 66-10273, 'Boron-Deoxidized Copper Withstands Brazing Temperatures' (1966)"
SMC045 = "Special Metals Corp., 'INCONEL alloy 718', Publication SMC-045 (2007)"
AK316 = "AK Steel, '316/316L Stainless Steel' product data bulletin"
SUMMERS = ("Summers, Chen, Rippe, Allen, Mouritz, Case, Lattimer, 'Overview of aluminum alloy mechanical "
           "properties during and after fires', Fire Science Reviews 4:3 (2015)")
EC9 = "EN 1999-1-2 (Eurocode 9), as reproduced in " + SUMMERS
EN1999 = "EN 1999-1-1 (Eurocode 9) section 3.2.5: nu = 0.3 for aluminium alloys"
CAO = ("Cao, Lin, Wang, Shi, Ma, Kang, Huang, 'Microstructure evolution and mechanical properties at high "
       "temperature of selective laser melted AlSi10Mg', J. Mater. Sci. Technol. 62 (2021) 162-172, Table 3")
EOS = "EOS GmbH, process data sheet 'EOS Aluminium AlSi10Mg for EOS M 290 | 30 um'"
MANSON = ("Manson, 'Fatigue: a complex subject - some simple approximations', Experimental Mechanics 5(7) "
          "(1965) 193-226: method of universal slopes")


# --------------------------------------------------------------------------------------
# Library
# --------------------------------------------------------------------------------------
def _cucrzr_e(t):   # UKAEA Eq. (1), average curve, SAA, T in °C, GPa
    return -1.9234e-4 * t * t - 2.1233e-2 * t + 124.91


def _cucrzr_ys(t):  # UKAEA Eq. (3), average curve, SAA, MPa
    return -2.2847e-4 * t * t - 0.13931 * t + 292.19


def _cucrzr_uts(t):  # UKAEA Eq. (5), average line, SAA, MPa
    return -0.42631 * t + 413.45


def _ec9_strain(t):  # EN 1999-1-2 thermal elongation of aluminium alloys, 20 <= T < 500 °C
    return 0.1e-7 * t * t + 22.5e-6 * t - 4.5e-4


_MATERIALS = [
    Material(
        name="CuCrZr", label="CuCrZr (copper alloy)", condition="solution annealed + aged (SAA), ITER grade",
        k=320.0, k_source=BAKER + " (293 K)",
        E=_tab(UKAEA + ", Eq. (1) average SAA fit, evaluated 20-500 °C (SAA data end near 500 °C)",
               [(t, _cucrzr_e(t)) for t in range(20, 501, 10)], GPA),
        alpha=_tab(BAKER, [(373.14, 16.3), (473.14, 17.0), (573.14, 17.6), (873.14, 18.3), (973.14, 19.0),
                           (1173.14, 19.8), (1356.14, 20.9)], 1e-6, "K"),
        nu=_tab(BAKER, [(293.14, 0.30), (873.14, 0.31)], 1.0, "K"),
        sy=_tab(UKAEA + ", Eq. (3) average SAA fit, evaluated 20-700 °C",
                [(t, _cucrzr_ys(t)) for t in range(20, 701, 10)], MPA),
        su=_tab(UKAEA + ", Eq. (5) average SAA fit, evaluated 20-700 °C",
                [(t, _cucrzr_uts(t)) for t in (20, 700)], MPA),
        lcf=LCF(2 * 3.5286e-2, -0.26199, (1e2, 1e6), 500 + C0,
                UKAEA + ", Eq. (7): Coffin-Manson fit, SAA, de/2 = 3.5286 Nf^-0.26199 %; data RT-500 °C")),
    Material(
        name="GRCop-42", label="GRCop-42 (copper alloy)",
        condition="wrought (extruded); Gradl et al. report L-PBF builds with strengths equivalent to wrought",
        k=346.0, k_source=GRADL + ", Fig. 1 (digitised from the vector figure, 23 °C)",
        E=_tab(MAHTABI + ": 123.5 GPa at room temperature (LP-DED, tensile modulus); held constant, no "
               "temperature data found", [(20, 123.5)], GPA),
        alpha=_secant(GRADL + ", Fig. 5 linear change (digitised), mean from 25 °C",
                      [(100, 0.1317e-2), (200, 0.3142e-2), (300, 0.4972e-2), (400, 0.6872e-2),
                       (500, 0.8682e-2), (600, 1.056e-2), (700, 1.256e-2), (800, 1.481e-2), (900, 1.726e-2),
                       (1000, 1.965e-2)], 25.0),
        nu=_tab("Pure copper value from " + NIST177 + " (295 K); no GRCop-specific data found", [(20, 0.34)]),
        sy=_tab(GRADL + ", Fig. 6 (digitised, wrought)",
                [(20, 199.9), (100, 186.1), (200, 171.0), (300, 156.3), (400, 140.6), (500, 121.8), (600, 98.4),
                 (700, 68.5), (800, 30.4)], MPA),
        su=_tab(GRADL + ", Fig. 7 (digitised, wrought)",
                [(20, 342.2), (100, 301.0), (200, 253.2), (300, 208.9), (400, 168.6), (500, 132.0), (600, 99.2),
                 (700, 70.3), (800, 45.1)], MPA),
        lcf=LCF(16.0772e-2, -0.33232, (561.0, 19402.0), 600 + C0,
                GRADL + ", Fig. 2: wrought GRCop-42 fit line, RT-600 °C pooled (digitised; axis labelled mm/mm, "
                "values are %)")),
    Material(
        name="Copper C101", label="Pure copper (C101)", condition="annealed oxygen-free copper",
        k=390.0, k_source="Typical handbook value kept from earlier versions of this tool (not re-verified)",
        E=_tab(NIST177 + ", Eq. (6-1) at 295 K; held constant above (data end at 300 K)", [(22, 126.0)], GPA),
        alpha=_tab(NIST177 + ", Eq. (7-4) at 293 K; held constant above (data end at 300 K)", [(20, 16.6)], 1e-6),
        nu=_tab(NIST177 + ", Eq. (6-6) at 295 K", [(22, 0.34)]),
        sy=_tab(NASA_TB + ": braze-coarsened OFHC, 5 ksi from room temperature to 900 °F",
                [(70, 5.0), (900, 5.0)], KSI, "F"),
        su=_tab(NASA_TB + ": 10 ksi at 900 °F, held constant at lower temperature (conservative)",
                [(900, 10.0)], KSI, "F")),
    Material(
        name="Al 6061-T6", label="Aluminium 6061", condition="T6 / T651",
        k=174.0, k_source=SUMMERS + ", Fig. 1, 6061-T651 at 25 °C (read from figure)",
        E=_tab(EC9 + ", Fig. 8 (read from figure)",
               [(20, 70), (100, 68), (150, 65), (200, 61), (250, 55), (300, 47.5), (350, 38), (400, 28),
                (450, 20)], GPA),
        alpha=_secant(EC9.replace("as reproduced in", "thermal elongation equation; checked against Fig. 3 of")
                      + ", mean from 20 °C",
                      [(t, _ec9_strain(t)) for t in range(30, 501, 10)], 20.0),
        nu=_tab(EN1999, [(20, 0.30)]),
        sy=_tab("Kaufman (2000) 6061-T651 data in " + SUMMERS + ", Fig. 9 (read from figure)",
                [(24, 275), (100, 260), (149, 250), (177, 240), (204, 220), (232, 195), (260, 150), (316, 75),
                 (371, 45), (427, 20)], MPA),
        su=_tab("Kaufman (2000) 6061-T651 data in " + SUMMERS + ", Fig. 10 (read from figure)",
                [(24, 310), (100, 285), (149, 260), (177, 250), (204, 235), (232, 200), (260, 160), (316, 88),
                 (371, 50), (427, 22)], MPA),
        ra=_tab(SUMMERS + ", Fig. 11, 6061-T651 reduction of area (read from figure)",
                [(25, 32), (100, 40), (150, 46), (200, 56), (250, 47), (300, 69), (350, 80), (400, 95)], 0.01)),
    Material(
        name="AlSi10Mg", label="AlSi10Mg (AM aluminium)", condition="laser powder bed fusion, as built",
        k=100.0, k_source=EOS + ": 100 W/m-K as manufactured, vertical (110 horizontal)",
        E=_tab(CAO, [(25, 76), (100, 65), (200, 59), (300, 49), (400, 18)], GPA),
        alpha=_tab(EOS + ": mean CTE 25-100, 25-200, 25-300 °C", [(100, 20.0), (200, 22.0), (300, 27.0)], 1e-6),
        nu=_tab(EN1999, [(20, 0.30)]),
        sy=_tab(CAO, [(25, 322), (100, 298), (200, 236), (300, 143), (400, 25)], MPA),
        su=_tab(CAO, [(25, 460), (100, 382), (200, 266), (300, 150), (400, 30)], MPA)),
    Material(
        name="Inconel 718", label="Inconel 718",
        condition="annealed 1800 °F + aged (strength: hot-rolled round, Table 20)",
        k=79 * BTU_IN, k_source=SMC045 + ", Table 5 (annealed + aged, 70 °F)",
        E=_tab(SMC045 + ", Table 4", [(70, 29.0), (200, 28.4), (400, 27.6), (600, 26.7), (800, 25.8),
                                      (1000, 24.8), (1200, 23.7), (1400, 22.3), (1600, 20.2)], KSI * 1e3, "F"),
        alpha=_tab(SMC045 + ", Table 5, mean from 70 °F",
                   [(200, 7.31), (400, 7.53), (600, 7.74), (800, 7.97), (1000, 8.09), (1200, 8.39), (1400, 8.91)],
                   1.8e-6, "F"),
        nu=_tab(SMC045 + ", Table 4", [(70, 0.294), (200, 0.288), (400, 0.280), (600, 0.272), (800, 0.271),
                                       (1000, 0.271), (1200, 0.283), (1400, 0.306), (1600, 0.331)], 1.0, "F"),
        sy=_tab(SMC045 + ", Table 20", [(80, 153.5), (1200, 145.0), (1300, 133.0)], KSI, "F"),
        su=_tab(SMC045 + ", Table 20", [(80, 190.5), (1200, 164.5), (1300, 145.5)], KSI, "F"),
        ra=_tab(SMC045 + ", Table 20, reduction of area", [(80, 32.5), (1200, 59.2), (1300, 34.0)], 0.01, "F")),
    Material(
        name="SS 316L", label="Stainless 316L", condition="annealed (bulletin covers 316 and 316L, typical values)",
        k=16.2, k_source=AK316 + " (value given at 100 °C)",
        E=_tab(AK316 + ": 193 GPa at room temperature; held constant, no temperature data in source",
               [(20, 193.0)], GPA),
        alpha=_tab(AK316 + ", mean from 0 °C", [(100, 16.0), (315, 16.2), (538, 17.5), (649, 18.5)], 1e-6),
        nu=_tab(AK316 + ": E/(2G) - 1 with E = 193 GPa, G = 77 GPa", [(20, 193.0 / (2 * 77.0) - 1)]),
        sy=_tab(AK316 + ", Tables 1-2", [(20, 290), (204, 241), (316, 214), (427, 190), (538, 165), (649, 145),
                                         (760, 124), (871, 110)], MPA),
        su=_tab(AK316 + ", Tables 1-2", [(20, 627), (204, 558), (316, 538), (427, 524), (538, 483), (649, 393),
                                         (760, 241), (871, 165)], MPA)),
]
MATERIALS: dict[str, Material] = {m.name: m for m in _MATERIALS}
BY_LABEL: dict[str, Material] = {m.label: m for m in _MATERIALS}
PROPS = {"E": "Young's modulus", "alpha": "thermal expansion", "nu": "Poisson's ratio",
         "sy": "yield strength", "su": "ultimate strength", "ra": "reduction of area"}


def get(name: str) -> Material:
    """Material by name or GUI label."""
    m = MATERIALS.get(name) or BY_LABEL.get(name)
    if m is None:
        raise ValueError(f"Unknown wall material '{name}'. Known: {', '.join(MATERIALS)}.")
    return m


def _warn(warnings: list[str] | None, msg: str) -> None:
    if warnings is not None and msg not in warnings:
        warnings.append(msg)


def prop(name: str, key: str, t: float, warnings: list[str] | None = None) -> float | None:
    """Property `key` (E, alpha, nu, sy, su) of a material at temperature t [K]; None if not tabulated."""
    m = get(name)
    tab: Table | None = getattr(m, key)
    if tab is None:
        return None
    if len(tab.t) == 1:
        _warn(warnings, f"{m.name}: {PROPS[key]} is a single value ({tab.t[0] - C0:.0f} °C), held constant "
                        "at all temperatures.")
        return tab.v[0]
    if t < tab.t[0] - 1e-6 or t > tab.t[-1] + 1e-6:
        _warn(warnings, f"{m.name}: {PROPS[key]} data cover {tab.t[0] - C0:.0f}-{tab.t[-1] - C0:.0f} °C; "
                        "the end value is used outside that range.")
    return float(np.interp(t, tab.t, tab.v))


def cycles_to_failure(name: str, d_eps: float, t: float, warnings: list[str] | None = None
                      ) -> tuple[float | None, str]:
    """Cycles to failure for a total strain range d_eps [-] at temperature t [K], and the method used.

    Measured curve if the material has one; otherwise Manson's universal slopes,
        d_eps = 3.5 (S_u / E) N^-0.12 + D^0.6 N^-0.6,   D = ln(1 / (1 - RA)),
    with S_u, E and RA at t (an estimate); None when neither is available.
    """
    m = get(name)
    if d_eps <= 0:
        return float("inf"), "no strain"
    if m.lcf is not None:
        c = m.lcf
        n = (d_eps / c.a) ** (1 / c.b)
        if not c.n_range[0] <= n <= c.n_range[1]:
            _warn(warnings, f"{m.name}: fatigue life extrapolated beyond the test data "
                            f"({c.n_range[0]:.0f}-{c.n_range[1]:.0f} cycles).")
        if t > c.t_max:
            _warn(warnings, f"{m.name}: fatigue data go up to {c.t_max - C0:.0f} °C; life at hotter points is "
                            "extrapolated.")
        return n, "measured"
    if m.ra is None or m.su is None:
        return None, "no fatigue data"
    e, su, ra = prop(name, "E", t, warnings), prop(name, "su", t, warnings), prop(name, "ra", t, warnings)
    d = np.log(1 / (1 - min(ra, 0.99)))
    f = lambda lg: 3.5 * su / e * 10 ** (-0.12 * lg) + d ** 0.6 * 10 ** (-0.6 * lg) - d_eps
    lo, hi = 0.0, 12.0
    if f(hi) > 0:
        return 10.0 ** hi, "universal slopes (estimate)"
    if f(lo) < 0:
        return 1.0, "universal slopes (estimate)"
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 10.0 ** (0.5 * (lo + hi)), "universal slopes (estimate)"


def summary_rows() -> list[tuple[str, str, float, str, str]]:
    """(name, condition, k, data range of yield strength, short sources) for reports and the README."""
    rows = []
    for m in _MATERIALS:
        rng = f"{m.sy.t[0] - C0:.0f}-{m.sy.t[-1] - C0:.0f} °C" if m.sy and len(m.sy.t) > 1 else "single value"
        rows.append((m.name, m.condition, m.k, rng, m.sy.source if m.sy else "-"))
    return rows
