"""
Isopropyl alcohol (2-propanol, IPA) coolant properties. CoolProp has no IPA model, so this module
builds one from published pure-component correlations, with the same interface as
regen_thermal.Coolant (name, p_crit, t_min, t_max, props, mu, h, t_from_h, t_sat, saturation).

Liquid (saturated-liquid correlations; the pressure effect on the compressed liquid is neglected)
    density            VDI Heat Atlas PPDS:  rho = rho_c + A tau^0.35 + B tau^(2/3) + C tau + D tau^(4/3)
    viscosity          VDI Heat Atlas PPDS eq. 9:  mu = E exp[A z^(1/3) + B z^(4/3)],  z = (C - T)/(T - D)
    thermal cond.      VDI Heat Atlas PPDS polynomial:  k = A + B T + C T^2 + D T^3 + E T^4
    heat capacity      Zabransky et al. (1996) recommended isobaric cubic splines, 185-473.2 K:
                       c_p / R = sum_j A_j+1 (T / 100 K)^j;  enthalpy is its analytical integral
Saturation
    vapour pressure    VDI Heat Atlas PPDS (Wagner 2.5-5 form):
                       ln(p/p_c) = (A tau + B tau^1.5 + C tau^2.5 + D tau^5) T_c / T
    latent heat        VDI Heat Atlas PPDS eq. 12:  h_fg = R T_c (A tau^(1/3) + B tau^(2/3) + C tau + D tau^2 + E tau^6)
Vapour (film model only)
    heat capacity      ideal gas, Poling, Prausnitz & O'Connell (2001) Appendix A polynomial
    thermal cond.      VDI Heat Atlas PPDS low-pressure gas polynomial
tau = 1 - T/T_c. Coefficients were transcribed from the data tables distributed with the `chemicals`
Python package (v1.5.2, C. Bell et al.) for CAS 67-63-0 and checked against handbook values (README).
"""

from __future__ import annotations

import math

R = 8.314462618        # J/mol-K
MW = 60.096            # g/mol (IUPAC critical-property compilation, as listed in `chemicals`)
TC, PC = 508.25, 4.762e6       # K, Pa (VDI Heat Atlas PPDS constants used with these fits)
T_MELT = 185.25                # K

RHO = (273.0, 865.9165, -744.0943, 975.9052, -381.1278)               # rho_c, A, B, C, D  (kg/m3)
MU_L = (3.07172, 3.84587, 580.414, -62.691, 7.42e-6)                  # A, B, C [K], D [K], E [Pa s]
K_L = (0.2028, -2.238e-4, -3.1e-8, 8.3e-11, -6.9e-14)                 # W/m-K, T in K
PSAT = (-8.44737, 1.17402, -6.97876, 0.69248)                          # Wagner 2.5-5
HVAP = (14.13231, -18.02015, 33.713095, -20.808299, -7.959708)         # PPDS 12
K_V = (0.001447, 1.1e-5, 1.2455e-7)                                    # W/m-K, T in K
CP_IG = (3.334, 0.018853, 3.644e-5, -6.115e-8, 2.543e-11)              # c_p/R, T in K
# Zabransky isobaric splines: (T_low, T_high, A1..A4)
CP_L = ((185.0, 270.0, (-11.2427, 32.7459, -15.5078, 2.65121)),
        (270.0, 360.0, (112.266, -104.486, 35.319, -3.6237)),
        (360.0, 473.2, (-175.513, 135.33, -31.2966, 2.5444)))
T_CP_DATA = CP_L[-1][1]        # K, top of the recommended c_p data
T_MAX = 0.98 * TC              # K, c_p extrapolated (with a warning) from T_CP_DATA to here
T_REF = 298.15                 # K, enthalpy zero

SOURCES = {
    "rho, mu, k (liquid), p_sat, h_fg, k (vapour)": "VDI-GVC (ed.), VDI Heat Atlas, 2nd ed., Springer 2010, "
                                                    "Section D3 (PPDS correlations)",
    "c_p (liquid)": "Zabransky, Ruzicka, Majer, Domalski, Heat Capacity of Liquids: Critical Review and "
                    "Recommended Values, J. Phys. Chem. Ref. Data Monograph 6, AIP 1996",
    "c_p (vapour, ideal gas)": "Poling, Prausnitz, O'Connell, The Properties of Gases and Liquids, 5th ed., "
                               "McGraw-Hill 2001, Appendix A",
    "transcription": "Coefficient tables of the `chemicals` Python package v1.5.2 (C. Bell et al.), CAS 67-63-0",
}


def _tau(t: float) -> float:
    return max(1.0 - t / TC, 0.0)


def density(t: float) -> float:
    rc, a, b, c, d = RHO
    x = _tau(t)
    return rc + a * x ** 0.35 + b * x ** (2 / 3) + c * x + d * x ** (4 / 3)


def viscosity(t: float) -> float:
    a, b, c, d, e = MU_L
    z = (c - t) / (t - d)
    z13 = math.copysign(abs(z) ** (1 / 3), z)
    return e * math.exp(a * z13 + b * z * z13)


def conductivity(t: float) -> float:
    return sum(c * t ** i for i, c in enumerate(K_L))


def _cp_piece(t: float):
    for lo, hi, a in CP_L:
        if t <= hi:
            return lo, a
    return CP_L[-1][0], CP_L[-1][2]


def cp_liquid(t: float) -> float:
    """J/kg-K"""
    _, a = _cp_piece(t)
    s = t / 100.0
    return R * (a[0] + a[1] * s + a[2] * s * s + a[3] * s ** 3) / (MW * 1e-3)


def _h_abs(t: float) -> float:
    """Integral of c_p from the first spline knot, J/kg (each piece integrated analytically; the last
    piece is extrapolated above T_CP_DATA)."""
    total = 0.0
    for lo, hi, a in CP_L:
        top = t if hi == T_CP_DATA else min(t, hi)
        if top > lo:
            f = lambda tt: 100.0 * sum(a[j] * (tt / 100.0) ** (j + 1) / (j + 1) for j in range(4))
            total += f(top) - f(lo)
        if t <= hi:
            break
    return R * total / (MW * 1e-3)


def enthalpy(t: float) -> float:
    """Liquid specific enthalpy relative to 298.15 K, J/kg."""
    return _h_abs(t) - _h_abs(T_REF)


def p_sat(t: float) -> float:
    a, b, c, d = PSAT
    x = _tau(t)
    return PC * math.exp((a * x + b * x ** 1.5 + c * x ** 2.5 + d * x ** 5) * TC / t)


def h_fg(t: float) -> float:
    """J/kg"""
    a, b, c, d, e = HVAP
    x = _tau(t)
    return R * TC * (a * x ** (1 / 3) + b * x ** (2 / 3) + c * x + d * x * x + e * x ** 6) / (MW * 1e-3)


def cp_vapour(t: float) -> float:
    """Ideal-gas c_p, J/kg-K."""
    return R * sum(c * t ** i for i, c in enumerate(CP_IG)) / (MW * 1e-3)


def k_vapour(t: float) -> float:
    return sum(c * t ** i for i, c in enumerate(K_V))


def _bisect(f, lo: float, hi: float, n: int = 80) -> float:
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


def t_saturation(p: float) -> float | None:
    if p >= PC:
        return None
    return _bisect(lambda t: p_sat(t) - p, 200.0, TC - 1e-6)


class IsopropanolCoolant:
    """IPA liquid with the regen_thermal.Coolant interface (correlations above; clamps to their range)."""

    def __init__(self, warnings: list[str]):
        self.name = "Isopropanol"
        self.warnings = warnings
        self.p_crit = PC
        self.t_min, self.t_max = T_MELT, T_MAX
        self._h_max = enthalpy(self.t_max)
        _warn(warnings, "Isopropanol coolant: built-in correlations (VDI Heat Atlas, Zabransky), "
                        "saturated-liquid properties with the pressure effect neglected. See README appendix.")

    def _t(self, t: float) -> float:
        if T_CP_DATA < t <= self.t_max:
            _warn(self.warnings, f"Isopropanol properties requested above {T_CP_DATA:.0f} K (top of the "
                                 "recommended c_p data): extrapolated toward the critical point, uncertain.")
        if t > self.t_max or t < self.t_min:
            _warn(self.warnings, f"Coolant property lookup outside the isopropanol correlation range "
                                 f"({self.t_min:.0f}-{self.t_max:.0f} K); clamped. Results near those "
                                 "stations are less reliable.")
            return min(max(t, self.t_min), self.t_max)
        return t

    def props(self, t: float, p: float) -> dict:
        t = self._t(t)
        rho, mu, cp, k = density(t), viscosity(t), cp_liquid(t), conductivity(t)
        return dict(rho=rho, mu=mu, cp=cp, k=k, pr=mu * cp / k)

    def mu(self, t: float, p: float) -> float:
        return viscosity(self._t(t))

    def h(self, t: float, p: float) -> float:
        return enthalpy(self._t(t))

    def t_from_h(self, h: float, p: float) -> float:
        if h >= self._h_max:
            _warn(self.warnings, "Coolant enthalpy exceeds the isopropanol correlation range: bulk "
                                 "temperature capped.")
            return self.t_max
        return _bisect(lambda t: enthalpy(t) - h, self.t_min, self.t_max)

    def t_sat(self, p: float) -> float | None:
        return t_saturation(p)

    def h_liquid(self, t: float, p: float) -> float:
        return enthalpy(self._t(t))

    def saturation(self, p: float) -> dict:
        """Film-model data at pressure p (see regen_thermal.Coolant.saturation)."""
        supercritical = p >= PC
        if supercritical:
            # No pseudo-critical model here: end the liquid region at the top of the property range
            t_end, hfg = self.t_max, 0.0
            _warn(self.warnings, f"Isopropanol above its critical pressure: the liquid film ends at "
                                 f"{t_end:.0f} K (0.98 T_c), the top of the property range.")
        else:
            t_end = min(t_saturation(p), self.t_max)
            hfg = h_fg(t_end)
        h_end = self.h_liquid(t_end, p)
        return dict(supercritical=supercritical, t_end=t_end, h_fg=hfg, h_end=h_end,
                    cp_v=cp_vapour(t_end), k_v=k_vapour(t_end), mw=MW * 1e-3)


def _warn(warnings: list[str], msg: str) -> None:
    if msg not in warnings:
        warnings.append(msg)
