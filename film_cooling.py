"""
Fuel film cooling: a non-entrained, non-reacting estimate, marched from the injector face.

Conventions (non-standard; see README):
    O/F is the CORE O/F (core oxidizer / core fuel).
    The film flow is a percentage of the CORE mass flow:  mdot_film = (pct / 100) * (mdot_ox + mdot_fuel,core).
    A 1 kg/s core with 0.15 kg/s of film is 15 %.

Model (Stechman-type two-region estimate):
1. Liquid film (Grisson, AEDC-TR-91-1, 1991, as reviewed by Giambelli et al., EUCASS 2023-648, sec. 2.2):
   non-reacting coolant, adiabatic film-wall interface, no entrainment mass loss. The film takes heat
   from the gas by convection only (no radiation, as elsewhere in this tool):
       heating:      d(T_L)/dA = h_g (T_aw - T_L) / (mdot_L c_p,L)            until T_L = T_sat(p)
       evaporation:  d(mdot_L)/dA = -h_b (T_aw - T_sat) / h_fg
   h_g is the local Bartz coefficient (with sigma at the film temperature). During evaporation the
   mass blowing reduces it, h_b = h_g ln(1 + B)/B with B = c_p,g (T_aw - T_sat)/h_fg (MW_g/MW_c)^0.6
   (transpiration correction; Giambelli et al. Eqs. 22, 31). The liquid film length L_c is where mdot_L = 0.
2. Gaseous film (Hatch & Papell, NASA TN D-130, 1959, Eq. 12; "no mixing", "absence of chemical
   reactions"), coolant = the evaporated fuel at T_sat from L_c on:
       ln(eta) = -[X - 0.04] (S V_g / alpha_c)^0.125 f(V_g/V_c),   X = integral(h_g dA) / (mdot_c c_p,c)
       eta = (T_aw - T_eff) / (T_aw - T_sat)
   with the vapour layer moving at the gas velocity (V_g/V_c = 1, f = 1, Hatch & Papell's optimum), so
   S V_g / alpha_c = mdot_c c_p,c / (2 pi r k_c) from continuity.
   Lucas & Golladay (NASA TN D-1988, 1963) found this correlation directly applicable in a small
   JP-4/GOX rocket chamber for a non-reactive coolant with convective heating; it needed a somewhat
   higher h in the convergent section and throat (film_h_scale anchors that).
3. Coupling to the wall (Stechman, Oberstone & Howell, J. Spacecraft Rockets 6(2), 1969, Eq. 1):
   q'' = h_g (T_eff - T_wg), i.e. the regen circuit keeps Bartz h_g and replaces the gas recovery
   temperature with the film driving temperature: T_L over the liquid film, T_eff downstream.

Liquid properties, saturation and latent heat are taken at chamber pressure. Supercritical films
(chamber pressure above the coolant's critical pressure) have no latent heat: the liquid region then
ends where the film reaches the pseudo-critical temperature, with a warning.
"""

from __future__ import annotations

import math

import numpy as np
from CoolProp.CoolProp import PropsSI

NONE, LIQUID, GAS = 0, 1, 2
REGION_NAMES = {NONE: "none", LIQUID: "liquid film", GAS: "gaseous film"}


def _warn(warnings: list[str], msg: str) -> None:
    if msg not in warnings:
        warnings.append(msg)


def velocity_ratio_factor(beta: float) -> float:
    """Hatch & Papell f(V_g/V_c): Eq. (10) for beta >= 1, Eq. (11) for beta <= 1 (angle in radians)."""
    if beta >= 1.0:
        return 1.0 + 0.4 * math.atan(beta - 1.0)
    r = 1.0 / beta
    return r ** (1.5 * (r - 1.0))


def hatch_papell_eta(x_param: float, sv_alpha: float, beta: float = 1.0) -> float:
    """Film effectiveness, Hatch & Papell Eq. (12); eta = 1 until x_param reaches 0.04."""
    if x_param <= 0.04:
        return 1.0
    return math.exp(-(x_param - 0.04) * sv_alpha ** 0.125 * velocity_ratio_factor(beta))


def pseudo_critical_temperature(fluid: str, p: float) -> float:
    """Temperature of maximum c_p at a supercritical pressure (CoolProp scan)."""
    t_c = PropsSI("Tcrit", fluid)
    ts = np.linspace(0.9 * t_c, 1.4 * t_c, 120)
    cps = []
    for t in ts:
        try:
            cps.append(PropsSI("C", "T", t, "P", p, fluid))
        except ValueError:
            cps.append(np.nan)
    return float(ts[int(np.nanargmax(cps))])


def film_profile(design, th, gas: dict, st: dict, cool, mdot_film: float, t_inj: float,
                 warnings: list[str]) -> dict:
    """March the film from the injector face. Returns per-station arrays and summary values.

    region (NONE / LIQUID / GAS), t_film (liquid temperature), mdot_liq, eta (gaseous region),
    t_drive (the temperature the wall sees: T_L or T_eff; NaN where there is no film), q_film
    (heat absorbed by the liquid at each station, W), and x_lc, l_c, q_liquid, t_end, supercritical.
    """
    from regen_thermal import bartz_h, recovery_temperature

    if not hasattr(cool, "p_crit"):
        raise ValueError("Film cooling needs a CoolProp coolant model (not a property table): saturation "
                         "and latent heat are taken from CoolProp.")
    fluid = cool.name
    x, r, ar, mach, da = st["x"], st["r"], st["ar"], st["mach"], st["da"]
    n = len(x)
    out = dict(region=np.full(n, NONE), t_film=np.full(n, np.nan), mdot_liq=np.full(n, np.nan),
               eta=np.full(n, np.nan), t_drive=np.full(n, np.nan), q_film=np.zeros(n))

    pc = gas["pc"]                                     # liquid properties and saturation at chamber pressure
    mw_ratio = (design.perf["mw_c"] / (PropsSI("M", fluid) * 1e3)) ** 0.6
    m_c = th.film_efficiency * mdot_film               # film that stays on the wall
    h_scale = th.hg_scale * th.film_h_scale

    supercritical = bool(pc >= cool.p_crit)
    if supercritical:
        _warn(warnings, f"Chamber pressure is above the critical pressure of {fluid} "
                        f"({cool.p_crit / 1e5:.1f} bar): the film is supercritical (no latent heat). Its liquid "
                        "region is taken to end at the pseudo-critical temperature; treat the film length as an "
                        "extrapolation of the model.")
        t_end = pseudo_critical_temperature(fluid, pc)
        h_fg = 0.0
    else:
        t_end = PropsSI("T", "P", pc, "Q", 0, fluid)
        h_fg = PropsSI("H", "P", pc, "Q", 1, fluid) - PropsSI("H", "P", pc, "Q", 0, fluid)
    if t_inj >= t_end:
        _warn(warnings, f"Film enters at {t_inj:.0f} K, at or above its saturation / pseudo-critical temperature "
                        f"({t_end:.0f} K): it flashes at injection.")

    t_l, m_l = min(t_inj, t_end), m_c
    # The liquid is marched in specific enthalpy (exact energy balance with variable c_p)
    h_end = PropsSI("H", "T", t_end - 0.01, "P", pc, fluid) if supercritical else PropsSI("H", "P", pc, "Q", 0, fluid)
    h_l = min(PropsSI("H", "T", t_l - 0.01, "P", pc, fluid), h_end)
    j_lc, x_lc, q_liquid = None, None, 0.0
    for j in range(n):
        t_aw = recovery_temperature(gas, mach[j])
        out["region"][j] = LIQUID
        out["t_film"][j], out["mdot_liq"][j], out["t_drive"][j] = t_l, m_l, t_l
        hg = h_scale * bartz_h(gas, ar[j], mach[j], t_l, th.throat_curv_ratio)
        q = max(hg * (t_aw - t_l) * da[j], 0.0)          # W of convection over this station's area
        need = m_l * (h_end - h_l)                       # heat to reach saturation
        if q < need:
            h_l += q / m_l
            t_l = PropsSI("T", "H", h_l, "P", pc, fluid)
            out["q_film"][j] = q
            q_liquid += q
            continue
        t_l, h_l = t_end, h_end
        if supercritical:
            out["q_film"][j] = need
            q_liquid += need
            j_lc, x_lc = j, x[j]
            break
        # Rest of the station evaporates the film; mass blowing lowers the convection (ln(1+B)/B)
        b = gas["cp"] * (t_aw - t_end) / h_fg * mw_ratio
        frac_blow = math.log(1 + b) / b if b > 1e-9 else 1.0
        q_evap = min((q - need) * frac_blow, m_l * h_fg)
        out["q_film"][j] = need + q_evap
        q_liquid += need + q_evap
        m_l -= q_evap / h_fg
        if m_l <= 1e-12 * m_c:
            j_lc, x_lc = j, x[j]
            break
    if j_lc is None:
        out["liquid_to_exit"] = True
    else:
        # Gaseous film downstream: the evaporated fuel at t_end, no mixing, no reaction
        out["region"][j_lc + 1:] = GAS
        if supercritical:
            cp_v = PropsSI("C", "T", t_end + 5.0, "P", pc, fluid)
            k_v = PropsSI("L", "T", t_end + 5.0, "P", pc, fluid)
        else:
            cp_v = PropsSI("C", "P", pc, "Q", 1, fluid)
            k_v = PropsSI("L", "P", pc, "Q", 1, fluid)
        sv_alpha = m_c * cp_v / (2 * math.pi * r[j_lc] * k_v)
        x_param = 0.0
        for j in range(j_lc + 1, n):
            t_aw = recovery_temperature(gas, mach[j])
            eta = hatch_papell_eta(x_param, sv_alpha)
            t_eff = t_aw - eta * (t_aw - t_end)
            out["eta"][j], out["t_drive"][j] = eta, t_eff
            hg = h_scale * bartz_h(gas, ar[j], mach[j], t_eff, th.throat_curv_ratio)
            x_param += hg * da[j] / (m_c * cp_v)
    out.update(x_lc=x_lc, l_c=None if x_lc is None else x_lc - x[0], q_liquid=q_liquid, t_end=t_end,
               supercritical=supercritical, mdot_film=mdot_film, mdot_on_wall=m_c, t_inj=t_inj,
               liquid_to_exit=j_lc is None)
    if j_lc is None:
        _warn(warnings, "The liquid film survives to the nozzle exit: the whole wall is film-protected.")
    return out
