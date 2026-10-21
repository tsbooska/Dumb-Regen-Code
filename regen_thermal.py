"""
Steady 1-D regenerative-cooling thermal circuit for the engine from shitty_rocket_code.py.

At each axial station:

    T_aw --[ 1/h_g (Bartz) ]--> T_wg --[ t_w/k_w ]--> T_wc --[ 1/h_c,eff (Nusselt + fins) ]--> T_b

    q''  = (T_aw - T_b) / (1/h_g + t_w/k_w + 1/h_c,eff)
    T_wg = T_aw - q''/h_g,     T_wc = T_wg - q'' t_w / k_w

Gas side
    Bartz:  h_g = 0.026/Dt^0.2 (mu^0.2 cp / Pr^0.6)_0 (Pc/c*)^0.8 (Dt/R_curv)^0.1 (At/A)^0.9 sigma
            sigma = 1 / ( [0.5 (T_wg/T0)(1 + (g-1)/2 M^2) + 0.5]^0.68 [1 + (g-1)/2 M^2]^0.12 )
    Recovery temperature  T_aw = T0 [1 + r (g-1)/2 M^2] / [1 + (g-1)/2 M^2],  r = Pr^(1/3)
    Frozen chamber transport properties from CEA. Local Mach from the area ratio.

Coolant side (rectangular channels, CoolProp properties at bulk T and P)
    Dh = 2wh/(w+h);  Sieder-Tate  Nu = 0.027 Re^0.8 Pr^(1/3) (mu_b/mu_w)^0.14
                  or Dittus-Boelter  Nu = 0.023 Re^0.8 Pr^0.4
    Ribs (lands) as straight fins:  land L = 2 pi (r + t_w)/N - w,  m = sqrt(2 h_c / (k_w L)),
        eta_f = tanh(m h)/(m h),   h_c,eff = h_c (w + 2 eta_f h) / (w + L)
    Energy balance marched in the coolant flow direction: dh_b = q'' dA / mdot_c
    Pressure drop: Darcy friction, Haaland friction factor.

Assumptions: steady, thin wall (planar conduction), constant k_wall, no radiation, no
boiling / film cooling, no structures. Stations outside the jacket are adiabatic (T_wg = T_aw).

Usage
-----
    py -3.14 regen_thermal.py                          # engine + thermal defaults from both dataclasses
    py -3.14 regen_thermal.py --pc 30 --of 2.3 --plot  # any shitty_rocket_code.py flag is accepted
    py -3.14 regen_thermal.py --channels channels.csv  # CSV columns: x_mm, width_mm, height_mm, wall_t_mm
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass, field

import numpy as np
from CoolProp.CoolProp import PropsSI

from shitty_rocket_code import (
    BAR, EngineDesign, _subsonic_mach, _supersonic_mach, design_engine, make_cea, parse_args,
    print_report,
)


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------
@dataclass
class ThermalInputs:
    # Coolant = the engine's fuel, at the full fuel flow, on its way to the injector
    coolant: str = "auto"                # "auto" = CoolProp fluid mapped from the CEA fuel name, or a CoolProp name
    t_in_K: float = 300.0                # fuel temperature entering the regen jacket
    p_in_bar: float = 40.0               # regen inlet pressure (injector.py solves this from the injector side)
    counterflow: bool = True             # True: enters at nozzle exit, exits at injector
    coolant_table: str | None = None     # CSV (T_K, rho, mu, cp, k) overriding the CoolProp fluid

    # Channels: rows of (x_mm, width_mm, height_mm, wall_t_mm), linearly interpolated in x
    # (x = 0 at the throat). One row = constant over the whole engine. With several rows the
    # jacket spans the table's x-range only, so end rows may sit past the injector / exit.
    # Station x positions come from the engine report and move if the engine is resized.
    n_channels: int = 40
    channel_table: list[tuple[float, float, float, float]] = field(
        default_factory=lambda: [
            #  x_mm,  width, height, wall_t   (mm)
            (-300.0,   1.5,   2.5,   1.0),   # cylinder (injector end; past the face is fine)
            ( -25.4,   1.5,   2.5,   1.0),   # cylinder end / start of convergence
            (   0.0,   0.9,   1.5,   0.7),   # throat
            (  50.0,   1.5,   2.5,   1.0),   # nozzle exit (past the exit is fine)
        ])

    # Wall
    k_wall: float = 150               # W/m-K (CuCrZr ~320, Inconel 718 ~15-25, SS 304 ~16-25)
    roughness_um: float = 6.0            # channel surface roughness
    t_wall_limit_K: float = 800.0        # hot-wall temperature limit, for warnings only

    # Gas side
    throat_curv_ratio: float = 1.0       # R_curv/Rt in Bartz (contour throat is sharp, so this is nominal)
    use_cstar_eff_T0: bool = True        # T0 = Tc * eta_cstar^2 instead of ideal Tc

    # Anchoring multipliers (e.g. 0.7 x Bartz from test data); 1.0 = correlation as published
    hg_scale: float = 1.0                # multiplies the Bartz gas-side coefficient h_g
    hc_scale: float = 1.0                # multiplies the coolant-side coefficient h_c (before fins)

    # Fuel film (amount set by EngineInputs.film_pct_core, % of CORE flow); see film_cooling.py
    film_efficiency: float = 1.0         # fraction of the film that stays on the wall (1 = no loss)
    film_h_scale: float = 1.0            # multiplies the gas-to-film h (Lucas & Golladay: higher near the throat)

    # Numerics / output
    correlation: str = "sieder-tate"     # "sieder-tate" or "dittus-boelter"
    n_stations: int = 300
    csv_path: str | None = "thermal_profile.csv"
    plot: bool = False
    plot_path: str | None = "thermal_profile.png"


# --------------------------------------------------------------------------------------
# Coolant properties
# --------------------------------------------------------------------------------------
# Fluids CoolProp has no model for, mapped to the usual surrogate
COOLANT_ALIASES = {"RP-1": "n-Dodecane", "RP1": "n-Dodecane", "Kerosene": "n-Dodecane",
                   "Jet-A": "n-Dodecane", "JetA": "n-Dodecane"}

# CEA fuel name -> CoolProp fluid used for coolant = "auto"
FUEL_COOLPROP = {"RP-1": "n-Dodecane", "RP1": "n-Dodecane", "JetA": "n-Dodecane", "Jet-A": "n-Dodecane",
                 "Ethanol": "Ethanol", "C2H5OH(L)": "Ethanol", "Methanol": "Methanol", "CH3OH(L)": "Methanol",
                 "CH4": "Methane", "LCH4_NASA": "Methane", "CH4(L)": "Methane", "LH2": "Hydrogen",
                 "GH2": "Hydrogen", "H2": "Hydrogen", "Propane": "Propane", "C3H8(L)": "Propane",
                 "NH3": "Ammonia", "Isopropanol": "IsopropylAlcohol"}


def make_coolant(design: EngineDesign, th: ThermalInputs, warnings: list[str]):
    """Property model for the fuel-as-coolant: table > explicit CoolProp name > auto from the fuel."""
    if th.coolant_table:
        return TableCoolant(th.coolant_table, warnings)
    name = th.coolant
    if name.lower() == "auto":
        fuel = design.inputs.fuel
        if fuel not in FUEL_COOLPROP:
            raise ValueError(f"No CoolProp model for fuel '{fuel}'. Supply a property table "
                             "(T_K, rho, mu, cp, k) or set the coolant to a CoolProp fluid name.")
        name = FUEL_COOLPROP[fuel]
        if name == "n-Dodecane":
            _warn(warnings, f"{fuel} is modeled with the CoolProp n-Dodecane surrogate. "
                            "Supply a property table for real data.")
    return Coolant(name, warnings)


class Coolant:
    """CoolProp wrapper that clamps temperatures to the fluid model's valid range."""

    def __init__(self, name: str, warnings: list[str]):
        if name in COOLANT_ALIASES:
            _warn(warnings, f"{name} is modeled with the CoolProp {COOLANT_ALIASES[name]} surrogate. "
                            "Supply a property table for real RP-1 data.")
            name = COOLANT_ALIASES[name]
        self.name = name
        self.warnings = warnings
        try:
            self.t_min = PropsSI("Tmin", name)
            self.t_max = PropsSI("Tmax", name)
            self.p_crit = PropsSI("pcrit", name)
        except ValueError as exc:
            raise ValueError(f"CoolProp does not know '{name}' ({exc}). For RP-1 use 'n-Dodecane'.") from None

    def _t(self, t: float) -> float:
        if t > self.t_max or t < self.t_min:
            _warn(self.warnings, f"Coolant property lookup outside CoolProp range for {self.name} "
                                 f"({self.t_min:.0f}-{self.t_max:.0f} K); clamped. Results near those "
                                 "stations are less reliable.")
            return min(max(t, self.t_min), self.t_max)
        return t

    def _get(self, q: str, t: float, p: float) -> float:
        try:
            return PropsSI(q, "T", t, "P", p, self.name)
        except ValueError:
            # Exactly on the saturation line (T, P don't fix the state): take the liquid side
            return PropsSI(q, "T", t - 0.05, "P", p, self.name)

    def props(self, t: float, p: float) -> dict:
        t = self._t(t)
        rho, mu, cp, k = (self._get(q, t, p) for q in ("D", "V", "C", "L"))
        return dict(rho=rho, mu=mu, cp=cp, k=k, pr=mu * cp / k)

    def mu(self, t: float, p: float) -> float:
        return self._get("V", self._t(t), p)

    def h(self, t: float, p: float) -> float:
        return self._get("H", self._t(t), p)

    def t_from_h(self, h: float, p: float) -> float:
        try:
            return PropsSI("T", "H", h, "P", p, self.name)
        except ValueError:
            _warn(self.warnings, "Coolant enthalpy exceeds the CoolProp range: bulk temperature capped.")
            return self.t_max

    def t_sat(self, p: float) -> float | None:
        return PropsSI("T", "P", p, "Q", 0, self.name) if p < self.p_crit else None


class TableCoolant:
    """Liquid properties vs temperature from a CSV (T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK).

    Pressure-independent; enthalpy is the integral of cp from the first table temperature.
    Same interface as Coolant.
    """

    def __init__(self, path: str, warnings: list[str]):
        tab = np.array(sorted(_read_numeric_csv(path, 5)), dtype=float)
        if len(tab) < 2:
            raise ValueError(f"Property table {path} needs at least 2 rows.")
        self.name = os.path.basename(path)
        self.warnings = warnings
        self.t, self.rho, self.mu_t, self.cp, self.k = tab.T
        self.t_min, self.t_max = self.t[0], self.t[-1]
        self.h_tab = np.concatenate([[0.0], np.cumsum(0.5 * (self.cp[1:] + self.cp[:-1]) * np.diff(self.t))])
        _warn(warnings, f"Coolant from table {self.name}: pressure-independent, boiling not checked.")

    def _t(self, t: float) -> float:
        if t > self.t_max or t < self.t_min:
            _warn(self.warnings, f"Coolant temperature outside the property table range "
                                 f"({self.t_min:.0f}-{self.t_max:.0f} K); clamped.")
            return min(max(t, self.t_min), self.t_max)
        return t

    def props(self, t: float, p: float) -> dict:
        t = self._t(t)
        rho, mu, cp, k = (float(np.interp(t, self.t, v)) for v in (self.rho, self.mu_t, self.cp, self.k))
        return dict(rho=rho, mu=mu, cp=cp, k=k, pr=mu * cp / k)

    def mu(self, t: float, p: float) -> float:
        return float(np.interp(self._t(t), self.t, self.mu_t))

    def h(self, t: float, p: float) -> float:
        return float(np.interp(self._t(t), self.t, self.h_tab))

    def t_from_h(self, h: float, p: float) -> float:
        if h > self.h_tab[-1]:
            _warn(self.warnings, "Coolant heated past the property table: bulk temperature capped.")
        return float(np.interp(h, self.h_tab, self.t))

    def t_sat(self, p: float) -> float | None:
        return None


def _warn(warnings: list[str], msg: str) -> None:
    if msg not in warnings:
        warnings.append(msg)


# --------------------------------------------------------------------------------------
# Gas side
# --------------------------------------------------------------------------------------
def gas_properties(design: EngineDesign, th: ThermalInputs) -> dict:
    """Frozen chamber transport properties and the reference values Bartz needs (SI)."""
    i, p = design.inputs, design.perf
    cea = make_cea(i.ox, i.fuel)
    cea_t = type(cea)(oxName=i.ox, fuelName=i.fuel, pressure_units="bar",
                      viscosity_units="poise", specific_heat_units="J/kg-K")
    cp, mu_poise, _, pr = (float(v) for v in
                           cea_t.get_Chamber_Transport(Pc=i.pc_bar, MR=p["of"], eps=p["eps"], frozen=1))
    t0 = p["tc"] * (i.eta_cstar ** 2 if th.use_cstar_eff_T0 else 1.0)
    return dict(cp=cp, mu=mu_poise * 0.1, pr=pr, t0=t0, gam=p["gam_t"],
                pc=i.pc_bar * BAR, cstar=p["cstar"], dt=design.sizing["dt"])


def bartz_h(gas: dict, area_ratio: float, mach: float, t_wg: float, curv_ratio: float) -> float:
    g, t0 = gas["gam"], gas["t0"]
    f = 1 + 0.5 * (g - 1) * mach ** 2
    sigma = 1.0 / ((0.5 * t_wg / t0 * f + 0.5) ** 0.68 * f ** 0.12)
    dt = gas["dt"]
    r_curv = curv_ratio * dt / 2
    return (0.026 / dt ** 0.2 * (gas["mu"] ** 0.2 * gas["cp"] / gas["pr"] ** 0.6)
            * (gas["pc"] / gas["cstar"]) ** 0.8 * (dt / r_curv) ** 0.1
            * area_ratio ** -0.9 * sigma)


def recovery_temperature(gas: dict, mach: float) -> float:
    g, r = gas["gam"], gas["pr"] ** (1 / 3)
    f = 1 + 0.5 * (g - 1) * mach ** 2
    return gas["t0"] * (1 + r * 0.5 * (g - 1) * mach ** 2) / f


# --------------------------------------------------------------------------------------
# Stations and channel geometry
# --------------------------------------------------------------------------------------
def build_stations(design: EngineDesign, n: int) -> dict:
    cx, cr = design.contour
    x = np.unique(np.concatenate([np.linspace(cx[0], cx[-1], n), cx]))
    r = np.interp(x, cx, cr)
    rt = design.sizing["rt"]
    ar = (r / rt) ** 2
    mach = np.array([_subsonic_mach(a, design.perf["gam_t"]) if xi < 0 else
                     (_supersonic_mach(a, design.perf["gam_t"]) if xi > 0 else 1.0)
                     for xi, a in zip(x, ar)])

    # Node-centred wetted area and slant length (half of each adjacent frustum segment)
    seg_s = np.hypot(np.diff(x), np.diff(r))
    seg_a = np.pi * (r[:-1] + r[1:]) * seg_s
    da = np.zeros_like(x)
    ds = np.zeros_like(x)
    da[:-1] += 0.5 * seg_a
    da[1:] += 0.5 * seg_a
    ds[:-1] += 0.5 * seg_s
    ds[1:] += 0.5 * seg_s
    return dict(x=x, r=r, ar=ar, mach=mach, da=da, ds=ds)


def channel_geometry(th: ThermalInputs, x: np.ndarray) -> tuple[np.ndarray, ...]:
    """Interpolate the channel table onto stations. Returns (cooled mask, w, h, t) in metres."""
    tab = np.array(sorted(th.channel_table), dtype=float)
    if len(tab) == 1:
        cooled = np.ones_like(x, dtype=bool)
    else:
        cooled = (x >= tab[0, 0] * 1e-3 - 1e-9) & (x <= tab[-1, 0] * 1e-3 + 1e-9)
    xs = tab[:, 0] * 1e-3
    w, h, t = (np.interp(x, xs, tab[:, k] * 1e-3) for k in (1, 2, 3))
    return cooled, w, h, t


def channel_table_from_sections(design: EngineDesign, cylinder, throat, exit_) -> list[tuple]:
    """Channel table from (width, height, wall_t) mm at three sections: constant along the
    cylinder, linear from the start of convergence to the throat, linear from throat to exit."""
    x_inj = design.chamber["x_inj"] * 1e3
    x_conv = -design.chamber["l_conv"] * 1e3
    x_exit = design.nozzle["length"] * 1e3
    return [(x_inj, *cylinder), (x_conv, *cylinder), (0.0, *throat), (x_exit, *exit_)]


def widths_for_ribs(design: EngineDesign, n_channels: int, rib_mm: tuple[float, float, float],
                    wall_t_mm: tuple[float, float, float]) -> tuple[float, float, float]:
    """Channel widths (mm) at cylinder, throat and exit that leave the given rib (land) widths.

    w = 2 pi (r + t_w) / N - rib at each section. Radius, wall thickness and rib all vary linearly
    between sections, so with channel_table_from_sections the rib follows a straight line from
    one section's value to the next (and is exactly constant if all three are equal).
    """
    s = design.sizing
    widths = tuple(2 * math.pi * (r * 1e3 + t) / n_channels - rib
                   for r, t, rib in zip((s["rc"], s["rt"], s["re"]), wall_t_mm, rib_mm))
    for name, w, r, t, rib in zip(("cylinder", "throat", "exit"), widths, (s["rc"], s["rt"], s["re"]),
                                  wall_t_mm, rib_mm):
        if w <= 0:
            pitch = 2 * math.pi * (r * 1e3 + t) / n_channels
            raise ValueError(f"A {rib:g} mm rib leaves no room for the channel at the {name} "
                             f"(pitch {pitch:.2f} mm with {n_channels} channels). "
                             "Use fewer channels or a thinner rib.")
    return widths


def _read_numeric_csv(path: str, ncols: int) -> list[tuple[float, ...]]:
    with open(path, newline="") as f:
        rows = [row for row in csv.reader(f) if row and not row[0].strip().startswith("#")]
    try:
        float(rows[0][0])
    except ValueError:
        rows = rows[1:]  # header
    return [tuple(float(v) for v in row[:ncols]) for row in rows]


def read_channel_csv(path: str) -> list[tuple[float, float, float, float]]:
    return _read_numeric_csv(path, 4)


# --------------------------------------------------------------------------------------
# Thermal circuit
# --------------------------------------------------------------------------------------
def haaland_f(re: float, rel_rough: float) -> float:
    if re < 2300:
        return 64.0 / re
    return (-1.8 * math.log10((rel_rough / 3.7) ** 1.11 + 6.9 / re)) ** -2


def solve_station(th, gas, cool, st, j, w, h, t, t_b, p_b, mdot_c, warnings, t_film=None) -> dict:
    """Converge T_wg / T_wc at one station for a given coolant bulk state.

    t_film: film driving temperature at this station (film_cooling.py); when given it replaces the gas
    recovery temperature in the circuit (Stechman: q'' = h_g (T_eff - T_wg)). t_aw is still reported.
    """
    ar, mach, r = st["ar"][j], st["mach"][j], st["r"][j]
    t_aw_gas = recovery_temperature(gas, mach)
    t_aw = t_film if t_film is not None and np.isfinite(t_film) else t_aw_gas
    n = th.n_channels
    dh = 2 * w * h / (w + h)
    land = 2 * math.pi * (r + t) / n - w
    if land <= 0:
        raise ValueError(f"Channels do not fit at x = {st['x'][j] * 1e3:.1f} mm: "
                         f"{n} x {w * 1e3:.2f} mm channels exceed the circumference. "
                         "Reduce n_channels or width.")

    cp_b = cool.props(t_b, p_b)
    vel = mdot_c / (n * w * h * cp_b["rho"])
    re = cp_b["rho"] * vel * dh / cp_b["mu"]

    t_wg = 0.5 * (t_aw + t_b)
    t_wc = t_b + 0.25 * (t_aw - t_b)
    for _ in range(200):
        hg = th.hg_scale * bartz_h(gas, ar, mach, t_wg, th.throat_curv_ratio)
        if th.correlation == "dittus-boelter":
            nu = 0.023 * re ** 0.8 * cp_b["pr"] ** 0.4
        else:
            nu = 0.027 * re ** 0.8 * cp_b["pr"] ** (1 / 3) * (cp_b["mu"] / cool.mu(t_wc, p_b)) ** 0.14
        hc = th.hc_scale * nu * cp_b["k"] / dh
        m = math.sqrt(2 * hc / (th.k_wall * land))
        eta_f = math.tanh(m * h) / (m * h)
        hce = hc * (w + 2 * eta_f * h) / (w + land)

        q = (t_aw - t_b) / (1 / hg + t / th.k_wall + 1 / hce)
        t_wg_new = t_aw - q / hg
        t_wc_new = t_wg_new - q * t / th.k_wall
        done = abs(t_wg_new - t_wg) < 0.1 and abs(t_wc_new - t_wc) < 0.1
        t_wg += 0.5 * (t_wg_new - t_wg)
        t_wc += 0.5 * (t_wc_new - t_wc)
        if done:
            break
    else:
        _warn(warnings, "Wall-temperature iteration did not fully converge at some stations.")

    f = haaland_f(re, th.roughness_um * 1e-6 / dh)
    dp = f * st["ds"][j] / dh * 0.5 * cp_b["rho"] * vel ** 2
    return dict(t_aw=t_aw_gas, t_hot=t_aw, hg=hg, hc=hc, eta_f=eta_f, q=q, t_wg=t_wg, t_wc=t_wc,
                vel=vel, re=re, dp=dp, land=land, dh=dh)


def run_thermal(design: EngineDesign, th: ThermalInputs, t_film_in: float | None = None) -> dict:
    """Regen march. With film (design.sizing['mdot_film'] > 0) the jacket carries core fuel + film fuel,
    and the film, injected at t_film_in (default: the regen inlet temperature; solve_feed_system passes
    the regen outlet temperature), sets the temperature the hot wall sees (film_cooling.py)."""
    warnings: list[str] = []
    gas = gas_properties(design, th)
    cool = make_coolant(design, th, warnings)
    st = build_stations(design, th.n_stations)
    cooled, w, h, t = channel_geometry(th, st["x"])
    mdot_film = design.sizing.get("mdot_film", 0.0)
    mdot_c = design.sizing["mdot_f"] + mdot_film

    film, film_warnings = None, []
    if mdot_film > 0:
        from film_cooling import film_profile
        film = film_profile(design, th, gas, st, cool, mdot_film,
                            th.t_in_K if t_film_in is None else t_film_in, film_warnings)
    t_film = film["t_drive"] if film else np.full(len(st["x"]), np.nan)

    n = len(st["x"])
    out = {k: np.full(n, np.nan) for k in
           ("t_aw", "t_hot", "hg", "hc", "eta_f", "q", "t_wg", "t_wc", "t_b", "p_b", "vel", "re", "land", "dh")}

    order = range(n - 1, -1, -1) if th.counterflow else range(n)
    t_b, p_b = th.t_in_K, th.p_in_bar * BAR
    h_b = cool.h(t_b, p_b)
    for j in order:
        if not cooled[j]:
            out["t_aw"][j] = recovery_temperature(gas, st["mach"][j])
            out["t_hot"][j] = out["t_wg"][j] = t_film[j] if np.isfinite(t_film[j]) else out["t_aw"][j]
            out["q"][j] = 0.0
            continue
        s = solve_station(th, gas, cool, st, j, w[j], h[j], t[j], t_b, p_b, mdot_c, warnings, t_film[j])
        for k in ("t_aw", "t_hot", "hg", "hc", "eta_f", "q", "t_wg", "t_wc", "vel", "re", "land", "dh"):
            out[k][j] = s[k]
        out["t_b"][j], out["p_b"][j] = t_b, p_b

        t_sat = cool.t_sat(p_b)
        if t_sat is not None and s["t_wc"] > t_sat:
            _warn(warnings, f"Cold wall exceeds coolant saturation T ({t_sat:.0f} K at "
                            f"{p_b / BAR:.1f} bar): nucleate boiling likely (not modeled).")
        if t_sat is not None and t_b > t_sat:
            _warn(warnings, "Coolant bulk temperature exceeds saturation: bulk boiling, results invalid.")

        h_b += s["q"] * st["da"][j] / mdot_c
        p_b -= s["dp"]
        if p_b <= 0:
            raise ValueError("Coolant pressure dropped below zero: raise p_in or enlarge the channels.")
        t_b = cool.t_from_h(h_b, p_b)

    re_min = np.nanmin(out["re"])
    if re_min < 1e4:
        _warn(warnings, f"Coolant Re drops to {re_min:.0f} (< 10^4): turbulent Nusselt correlation is "
                        "outside its validity range there. Use smaller channels or more coolant flow.")
    if np.nanmax(out["t_wg"][cooled]) > th.t_wall_limit_K:
        _warn(warnings, f"Hot-wall temperature exceeds the {th.t_wall_limit_K:.0f} K limit.")
    _warn(warnings, f"Contour throat is sharp: the throat heat-flux peak depends on throat_curv_ratio "
                    f"= {th.throat_curv_ratio:g} (Bartz curvature term).")
    for msg in film_warnings:
        _warn(warnings, msg)

    q_total = float(np.sum(out["q"] * st["da"]))
    return dict(st=st, cooled=cooled, w=w, h=h, t=t, gas=gas, mdot_c=mdot_c, out=out, cool=cool,
                t_out=t_b, p_out=p_b, dp_regen=th.p_in_bar * BAR - p_b, q_total=q_total, warnings=warnings,
                film=film)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def print_thermal_report(design: EngineDesign, th: ThermalInputs, res: dict) -> None:
    st, o, gas, cooled = res["st"], res["out"], res["gas"], res["cooled"]
    x = st["x"]
    jmax = int(np.nanargmax(np.where(cooled, o["t_wg"], np.nan)))

    print("=" * 70)
    print(f" REGEN THERMAL: {design.inputs.fuel} ({res['cool'].name}) in {th.n_channels} channels, "
          f"{'counterflow' if th.counterflow else 'coflow'}, k_wall = {th.k_wall:g} W/m-K")
    print("=" * 70)
    print(f"  Gas (frozen, chamber): cp = {gas['cp']:.0f} J/kg-K, mu = {gas['mu']:.3e} Pa-s, "
          f"Pr = {gas['pr']:.3f}, T0 = {gas['t0']:.0f} K")
    if th.hg_scale != 1.0 or th.hc_scale != 1.0:
        print(f"  Anchoring multipliers    h_g x {th.hg_scale:g} (Bartz), h_c x {th.hc_scale:g} (coolant)")
    film = res.get("film")
    print(f"  Coolant flow             {res['mdot_c']:10.4f} kg/s  "
          f"({'= core fuel + film fuel' if film else '= engine fuel flow'})")
    print(f"  Fuel T regen in / out    {th.t_in_K:10.1f} / {res['t_out']:.1f} K   "
          f"(dT = {res['t_out'] - th.t_in_K:.1f} K)")
    print(f"  Fuel P regen in / out    {th.p_in_bar:10.2f} / {res['p_out'] / BAR:.2f} bar   "
          f"(dP = {res['dp_regen'] / BAR:.2f} bar)")
    print(f"  Total heat load          {res['q_total'] / 1e3:10.2f} kW")
    print(f"  Max hot-wall T           {o['t_wg'][jmax]:10.1f} K at x = {x[jmax] * 1e3:.1f} mm")
    print(f"  Max heat flux            {np.nanmax(o['q']) / 1e6:10.2f} MW/m^2")
    if film:
        print_film_summary(design, res)

    keys = [("injector", 0), ("conv. start", int(np.argmin(abs(x + design.chamber["l_conv"])))),
            ("throat", int(np.argmin(abs(x)))), ("exit", len(x) - 1), ("max T_wg", jmax)]
    print(f"\n  {'station':<12}{'x mm':>8}{'A/At':>7}{'M':>6}{'T_aw':>7}{'h_g':>7}{'h_c':>7}"
          f"{'eta_f':>6}{'q MW/m2':>8}{'T_wg':>7}{'T_wc':>7}{'T_b':>7}{'V m/s':>7}")
    print(f"  {'':<12}{'':>8}{'':>7}{'':>6}{'K':>7}{'kW/m2K':>7}{'kW/m2K':>7}{'':>6}{'':>8}"
          f"{'K':>7}{'K':>7}{'K':>7}{'':>7}")
    for name, j in keys:
        print(f"  {name:<12}{x[j] * 1e3:8.1f}{st['ar'][j]:7.2f}{st['mach'][j]:6.2f}{o['t_aw'][j]:7.0f}"
              f"{o['hg'][j] / 1e3:7.2f}{o['hc'][j] / 1e3:7.1f}{o['eta_f'][j]:6.2f}{o['q'][j] / 1e6:8.2f}"
              f"{o['t_wg'][j]:7.0f}{o['t_wc'][j]:7.0f}{o['t_b'][j]:7.0f}{o['vel'][j]:7.2f}")

    if res["warnings"]:
        print("\n WARNINGS")
        for msg in res["warnings"]:
            print(f"  ! {msg}")
    print()


def print_film_summary(design: EngineDesign, res: dict) -> None:
    """Film-cooling lines for the thermal report (only printed when a film is used)."""
    f, x = res["film"], res["st"]["x"]
    print(f"  Film coolant             {f['mdot_film']:10.4f} kg/s  ({design.inputs.film_pct_core:g} % of core flow; "
          f"{f['mdot_on_wall']:.4f} kg/s stays on the wall), injected at {f['t_inj']:.0f} K")
    if f["liquid_to_exit"]:
        print("  LIQUID FILM LENGTH       whole engine (the liquid film reaches the nozzle exit)")
    else:
        print(f"  LIQUID FILM LENGTH       {f['l_c'] * 1e3:10.1f} mm   from the injector face "
              f"(ends at x = {f['x_lc'] * 1e3:.1f} mm)")
    print(f"  Film {'pseudo-critical' if f['supercritical'] else 'saturation'} T       {f['t_end']:10.1f} K   "
          f"(at Pc; heat absorbed by the liquid film {f['q_liquid'] / 1e3:.2f} kW)")
    gas_j = np.flatnonzero(np.isfinite(f["eta"]))
    if len(gas_j):
        j_t = int(np.argmin(abs(x)))
        eta_t = f["eta"][j_t] if np.isfinite(f["eta"][j_t]) else 1.0
        print(f"  Gaseous film (Hatch-Papell) effectiveness: {f['eta'][gas_j[0]]:.2f} where it starts, "
              f"{eta_t:.2f} at the throat, {f['eta'][-1]:.2f} at the exit")


def write_thermal_csv(res: dict, path: str) -> None:
    st, o = res["st"], res["out"]
    cols = [("x_mm", st["x"] * 1e3), ("r_mm", st["r"] * 1e3), ("A_At", st["ar"]), ("Mach", st["mach"]),
            ("chan_w_mm", res["w"] * 1e3), ("chan_h_mm", res["h"] * 1e3), ("wall_t_mm", res["t"] * 1e3),
            ("land_mm", o["land"] * 1e3), ("Dh_mm", o["dh"] * 1e3),
            ("T_aw_K", o["t_aw"]), ("h_g_W_m2K", o["hg"]), ("h_c_W_m2K", o["hc"]), ("eta_fin", o["eta_f"]),
            ("q_W_m2", o["q"]), ("T_wg_K", o["t_wg"]), ("T_wc_K", o["t_wc"]), ("T_b_K", o["t_b"]),
            ("P_b_bar", o["p_b"] / BAR), ("V_c_m_s", o["vel"]), ("Re_c", o["re"])]
    if res.get("film"):
        f = res["film"]
        cols += [("film_region", f["region"].astype(float)), ("T_film_liquid_K", f["t_film"]),
                 ("eta_film_gas", f["eta"]), ("T_drive_film_K", o["t_hot"])]
    with open(path, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow([c for c, _ in cols])
        for row in zip(*(v for _, v in cols)):
            wr.writerow([f"{v:.6g}" for v in row])
    print(f"Thermal profile written to {path}  ({len(st['x'])} stations)")


def plot_thermal(design: EngineDesign, th: ThermalInputs, res: dict, show: bool) -> None:
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    st, o = res["st"], res["out"]
    x = st["x"] * 1e3
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                                     gridspec_kw=dict(height_ratios=[3, 2, 1.3]))
    a1.plot(x, o["t_wg"], "tab:red", label="T_wg (hot wall)")
    a1.plot(x, o["t_wc"], "tab:orange", label="T_wc (cold wall)")
    a1.plot(x, o["t_b"], "tab:blue", label="T_b (coolant)")
    a1.axhline(th.t_wall_limit_K, color="k", ls=":", lw=1, label="wall limit")
    a1.set_ylabel("T [K]")
    a1.legend(loc="upper left")
    a1.grid(alpha=0.3)

    a2.plot(x, o["q"] / 1e6, "tab:purple", label="q'' [MW/m²]")
    a2.set_ylabel("q'' [MW/m²]")
    a2b = a2.twinx()
    a2b.plot(x, o["hg"] / 1e3, "tab:green", ls="--", label="h_g [kW/m²K]")
    a2b.set_ylabel("h_g [kW/m²K]")
    a2.grid(alpha=0.3)
    a2.legend(loc="upper left")
    a2b.legend(loc="upper right")

    a3.plot(x, st["r"] * 1e3, "k")
    a3.fill_between(x, 0, st["r"] * 1e3, color="tab:orange", alpha=0.12)
    a3.set_ylabel("r [mm]")
    a3.set_xlabel("x [mm]")
    a3.grid(alpha=0.3)

    i = design.inputs
    a1.set_title(f"{i.ox}/{i.fuel} Pc={i.pc_bar:g} bar  F={i.thrust_n:g} N  |  {res['cool'].name}, "
                 f"{th.n_channels} channels, k_w={th.k_wall:g} W/m-K")
    fig.tight_layout()
    if th.plot_path:
        fig.savefig(th.plot_path, dpi=150)
        print(f"Plot written to {th.plot_path}")
    if show:
        plt.show()


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv=None) -> None:
    d = ThermalInputs()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)
    ap.add_argument("--channels", help="channel table CSV: x_mm, width_mm, height_mm, wall_t_mm")
    ap.add_argument("--coolant-table", default=d.coolant_table,
                    help="coolant property CSV: T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK")
    ap.add_argument("--thermal-csv", default=d.csv_path, help="thermal profile CSV path ('' to skip)")
    ap.add_argument("--thermal-plot-path", default=d.plot_path)
    o, engine_argv = ap.parse_known_args(argv)

    inp, show = parse_args(engine_argv)
    th = ThermalInputs(csv_path=o.thermal_csv or None, plot=inp.plot or d.plot,
                       plot_path=o.thermal_plot_path, coolant_table=o.coolant_table)
    if o.channels:
        th.channel_table = read_channel_csv(o.channels)

    design = design_engine(inp)
    print_report(design)
    res = run_thermal(design, th)
    print_thermal_report(design, th, res)
    if th.csv_path:
        write_thermal_csv(res, th.csv_path)
    if th.plot:
        plot_thermal(design, th, res, show)


if __name__ == "__main__":
    main()
