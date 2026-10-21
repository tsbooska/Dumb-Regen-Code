"""
Regen channel optimizer: number of channels and channel width / height at the cylinder, throat and
nozzle exit that minimise the peak hot-wall temperature.

    minimise    max T_wg (cooled stations)
    over        N (integer),  w_cyl, w_thr, w_exit,  h_cyl, h_thr, h_exit
    subject to  w >= w_min,   h_min <= h <= h_max,   rib (land) >= rib_min,   dP_regen <= dP_max

Wall thickness is held at the given section values. The land 2 pi (r + t_w)/N - w is linear between
sections (r, t_w and w all are), so the rib limit is enforced exactly as an upper bound on w at the
three sections:  w <= 2 pi (r + t_w)/N - rib_min.

Method
    Inner (fixed N): COBYQA (derivative-free trust region, scipy) on the six sizes scaled to [0, 1],
        with dP_regen <= dP_max as a nonlinear constraint. run_thermal is evaluated at
        p_in = Pc (1 + s_f) + dP_max with n_stations_opt stations.
    Outer (N): coarse scan n_min..n_max in steps of n_step, warm-starting each N from the previous
        optimum, then a step-1 refinement around the best N.
    The optimum is re-run at the full station count through solve_feed_system.

Usage
-----
    py -3.14 regen_optimize.py                                 # engine, thermal, injector defaults
    py -3.14 regen_optimize.py --dp-max 3 --rib-min 1.0        # any shitty_rocket_code.py flag too
"""

from __future__ import annotations

import argparse
import dataclasses
from dataclasses import dataclass

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize

import injector as inj
import regen_thermal as rt
from shitty_rocket_code import BAR, EngineDesign, design_engine, parse_args, print_report

SECTIONS = ("cylinder", "throat", "exit")
PENALTY_T, PENALTY_DP = 1.0e4, 1.0e3     # objective / constraint for a run_thermal that fails


# --------------------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------------------
@dataclass
class OptimizerInputs:
    n_min: int = 20                      # channel count range searched
    n_max: int = 80
    n_step: int = 4                      # coarse-scan step in N (refined with step 1 afterwards)
    w_min_mm: float = 0.8                # minimum channel width
    h_min_mm: float = 1.0                # minimum channel height
    h_max_mm: float = 6.0                # maximum channel height
    rib_min_mm: float = 0.8              # minimum rib (land) width
    dp_max_bar: float = 5.0              # allowed regen channel pressure drop
    n_stations_opt: int = 100            # axial stations used while optimising (final check uses th's)
    max_evals_per_n: int = 60            # run_thermal calls per channel count


def sections_from_table(th: rt.ThermalInputs) -> tuple[tuple[float, float, float], ...]:
    """(width, height, wall_t) mm at the cylinder, throat and exit rows of a channel table."""
    tab = sorted(th.channel_table)
    return tab[0][1:], min(tab, key=lambda row: abs(row[0]))[1:], tab[-1][1:]


# --------------------------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------------------------
def _with_channels(design: EngineDesign, th: rt.ThermalInputs, n: int, x_mm, wall_t_mm) -> rt.ThermalInputs:
    sections = [(x_mm[k], x_mm[3 + k], wall_t_mm[k]) for k in range(3)]
    return dataclasses.replace(th, n_channels=n, channel_table=rt.channel_table_from_sections(design, *sections))


def _peak_t_wg(res: dict) -> float:
    return float(np.nanmax(np.where(res["cooled"], res["out"]["t_wg"], np.nan)))


def _optimize_n(design, th, wall_t_mm, n, opt, lo, hi, x_start_mm, radius) -> dict:
    """COBYQA on the six sizes for one channel count. Returns the best feasible point it evaluated."""
    cache: dict[tuple, tuple[float, float]] = {}

    def evaluate(u):
        u = np.clip(u, 0.0, 1.0)
        key = tuple(np.round(u, 6))
        if key not in cache:
            x = lo + u * (hi - lo)
            try:
                res = rt.run_thermal(design, _with_channels(design, th, n, x, wall_t_mm))
                cache[key] = (_peak_t_wg(res), res["dp_regen"] / BAR)
            except ValueError:  # coolant pressure ran out, or the channels do not fit
                cache[key] = (PENALTY_T, PENALTY_DP)
        return cache[key]

    u0 = np.clip((np.asarray(x_start_mm) - lo) / (hi - lo), 0.0, 1.0)
    minimize(lambda u: evaluate(u)[0], u0, method="COBYQA", bounds=[(0.0, 1.0)] * 6,
             constraints=[NonlinearConstraint(lambda u: evaluate(u)[1], -np.inf, opt.dp_max_bar)],
             options=dict(maxfev=opt.max_evals_per_n, initial_tr_radius=radius, final_tr_radius=1e-3))

    feasible = [(t, dp, u) for u, (t, dp) in cache.items() if dp <= opt.dp_max_bar and t < PENALTY_T]
    if not feasible:
        return dict(n=n, feasible=False, evals=len(cache))
    t, dp, u = min(feasible)
    return dict(n=n, feasible=True, t_wg=t, dp_bar=dp, x_mm=lo + np.array(u) * (hi - lo), evals=len(cache))


# --------------------------------------------------------------------------------------
# Optimizer
# --------------------------------------------------------------------------------------
def optimize_channels(design: EngineDesign, th: rt.ThermalInputs, inj_inp: inj.InjectorInputs,
                      wall_t_mm: tuple[float, float, float], opt: OptimizerInputs,
                      start_mm: tuple[float, ...] | None = None, progress=None) -> dict:
    """Best N and (w, h) at the three sections. start_mm = (w_cyl, w_thr, w_exit, h_cyl, h_thr, h_exit)
    is the first starting point (default: mid-range)."""
    say = progress or (lambda msg: None)
    if opt.n_min < 3 or opt.n_max < opt.n_min or opt.n_step < 1:
        raise ValueError("Optimizer: need 3 <= N min <= N max and N step >= 1.")
    if opt.h_max_mm <= opt.h_min_mm:
        raise ValueError("Optimizer: the maximum channel height must exceed the minimum.")
    p_inj_f = design.inputs.pc_bar * BAR * (1 + inj_inp.stiff_fuel)
    th_opt = dataclasses.replace(th, n_stations=opt.n_stations_opt, p_in_bar=p_inj_f / BAR + opt.dp_max_bar,
                                 csv_path=None, plot=False)

    def bounds(n):
        try:
            w_max = np.array(rt.widths_for_ribs(design, n, (opt.rib_min_mm,) * 3, wall_t_mm))
        except ValueError:
            return None
        if (w_max <= opt.w_min_mm).any():
            return None
        return (np.concatenate([[opt.w_min_mm] * 3, [opt.h_min_mm] * 3]),
                np.concatenate([w_max, [opt.h_max_mm] * 3]))

    history: dict[int, dict] = {}
    x_warm = None if start_mm is None else np.asarray(start_mm, dtype=float)

    def run(n, x_start, radius):
        b = bounds(n)
        if b is None:
            history[n] = dict(n=n, feasible=False, evals=0, reason="rib limit")
            return None
        lo, hi = b
        if x_start is None:
            x_start = 0.5 * (lo + hi)
        r = _optimize_n(design, th_opt, wall_t_mm, n, opt, lo, hi, x_start, radius)
        history[n] = r
        if r["feasible"]:
            say(f"Optimizing channels: N = {n}, peak T_wg {r['t_wg']:.0f} K, dP {r['dp_bar']:.2f} bar")
        else:
            say(f"Optimizing channels: N = {n}, no layout meets the dP limit")
        return r

    # Coarse scan. The rib bound tightens as N grows, so the first N that cannot fit ends the scan.
    for n in range(opt.n_min, opt.n_max + 1, opt.n_step):
        r = run(n, x_warm, 0.25 if x_warm is None else 0.1)
        if r is None:
            break
        if r["feasible"]:
            x_warm = r["x_mm"]

    def best():
        ok = [r for r in history.values() if r["feasible"]]
        return min(ok, key=lambda r: r["t_wg"]) if ok else None

    b0 = best()
    if b0 is None:
        raise ValueError("No feasible channel layout: every N either violates the rib limit or exceeds the "
                         "regen dP limit. Raise dP max, lower the minimum sizes, or widen the N range.")
    # Step-1 refinement around the best coarse N
    for n in range(max(opt.n_min, b0["n"] - opt.n_step + 1), min(opt.n_max, b0["n"] + opt.n_step - 1) + 1):
        if n not in history:
            run(n, b0["x_mm"], 0.1)
    b = best()

    # Verify at the full station count with the secant-converged feed system
    th_best = _with_channels(design, th, b["n"], b["x_mm"], wall_t_mm)
    th_best = dataclasses.replace(th_best, csv_path=None, plot=False)
    feed = inj.solve_feed_system(design, th_best, inj_inp)
    res = feed["thermal"]
    t_wg, dp = _peak_t_wg(res), res["dp_regen"] / BAR

    warnings: list[str] = []
    x = b["x_mm"]
    if dp > opt.dp_max_bar * 1.001:
        warnings.append(f"Regen dP at {th.n_stations} stations is {dp:.2f} bar, above the {opt.dp_max_bar:g} bar "
                        f"limit met at {opt.n_stations_opt} stations. Raise the optimizer's station count.")
    if np.isclose(x[3:], opt.h_max_mm, rtol=0.01).any():
        warnings.append(f"Channel height reached its {opt.h_max_mm:g} mm upper bound; a larger h max may "
                        "lower T_wg further.")
    if b["n"] in (opt.n_min, opt.n_max):
        warnings.append(f"Best N = {b['n']} is at the edge of the searched range; widen N min / N max.")
    elif history.get(b["n"] + 1, {}).get("reason") == "rib limit":
        warnings.append(f"Best N = {b['n']} is the most channels the {opt.rib_min_mm:g} mm rib and "
                        f"{opt.w_min_mm:g} mm width limits allow: the manufacturing limits, not dP, set the optimum.")
    lands = res["out"]["land"][res["cooled"]] * 1e3
    if np.nanmin(lands) < opt.rib_min_mm - 1e-6:
        warnings.append(f"Minimum rib {np.nanmin(lands):.3f} mm is below the {opt.rib_min_mm:g} mm limit.")

    sections = tuple((x[k], x[3 + k], wall_t_mm[k]) for k in range(3))
    return dict(n=b["n"], sections=sections, t_wg=t_wg, dp_bar=dp, rib_min=float(np.nanmin(lands)),
                th=th_best, feed=feed, history=[history[n] for n in sorted(history)],
                evals=sum(r["evals"] for r in history.values()), warnings=warnings, inputs=opt)


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------
def print_optimize_report(opt_res: dict, baseline: dict | None = None) -> None:
    o = opt_res["inputs"]
    print("=" * 70)
    print(" CHANNEL OPTIMIZER: minimise peak hot-wall temperature")
    print("=" * 70)
    print(f"  Limits: w >= {o.w_min_mm:g} mm, {o.h_min_mm:g} <= h <= {o.h_max_mm:g} mm, "
          f"rib >= {o.rib_min_mm:g} mm, dP_regen <= {o.dp_max_bar:g} bar")
    print(f"\n  {'N':>4}{'peak T_wg K':>13}{'dP bar':>9}{'evals':>7}")
    for r in opt_res["history"]:
        if r["feasible"]:
            mark = "  <- best" if r["n"] == opt_res["n"] else ""
            print(f"  {r['n']:4d}{r['t_wg']:13.1f}{r['dp_bar']:9.2f}{r['evals']:7d}{mark}")
        else:
            print(f"  {r['n']:4d}{'infeasible':>13}{'':>9}{r['evals']:7d}  ({r.get('reason', 'dP limit')})")
    print(f"\n  BEST: {opt_res['n']} channels   ({opt_res['evals']} thermal runs)")
    print(f"  {'section':<10}{'width':>8}{'height':>8}{'wall t':>8}   mm")
    for name, (w, h, t) in zip(SECTIONS, opt_res["sections"]):
        print(f"  {name:<10}{w:8.3f}{h:8.3f}{t:8.3f}")
    print(f"  Peak T_wg {opt_res['t_wg']:.1f} K, regen dP {opt_res['dp_bar']:.2f} bar, "
          f"minimum rib {opt_res['rib_min']:.3f} mm (full station count)")
    if baseline:
        print(f"  Starting layout: peak T_wg {baseline['t_wg']:.1f} K, regen dP {baseline['dp_bar']:.2f} bar")
    if opt_res["warnings"]:
        print("\n WARNINGS")
        for msg in opt_res["warnings"]:
            print(f"  ! {msg}")
    print()


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv=None) -> None:
    d = OptimizerInputs()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                 add_help=False)
    a = ap.add_argument
    a("--channels", help="starting channel table CSV: x_mm, width_mm, height_mm, wall_t_mm")
    a("--coolant-table", help="fuel property CSV: T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK")
    a("--n-min", type=int, default=d.n_min, help="smallest channel count searched")
    a("--n-max", type=int, default=d.n_max, help="largest channel count searched")
    a("--n-step", type=int, default=d.n_step, help="coarse-scan step in N")
    a("--w-min", type=float, default=d.w_min_mm, help="minimum channel width [mm]")
    a("--h-min", type=float, default=d.h_min_mm, help="minimum channel height [mm]")
    a("--h-max", type=float, default=d.h_max_mm, help="maximum channel height [mm]")
    a("--rib-min", type=float, default=d.rib_min_mm, help="minimum rib (land) width [mm]")
    a("--dp-max", type=float, default=d.dp_max_bar, help="allowed regen channel pressure drop [bar]")
    a("--opt-stations", type=int, default=d.n_stations_opt, help="axial stations while optimising")
    a("--evals-per-n", type=int, default=d.max_evals_per_n, help="thermal runs per channel count")
    o, engine_argv = ap.parse_known_args(argv)

    inp, _ = parse_args(engine_argv)
    inp.csv_path = None
    th = rt.ThermalInputs(csv_path=None, plot=False, coolant_table=o.coolant_table)
    if o.channels:
        th.channel_table = rt.read_channel_csv(o.channels)
    opt = OptimizerInputs(n_min=o.n_min, n_max=o.n_max, n_step=o.n_step, w_min_mm=o.w_min, h_min_mm=o.h_min,
                          h_max_mm=o.h_max, rib_min_mm=o.rib_min, dp_max_bar=o.dp_max,
                          n_stations_opt=o.opt_stations, max_evals_per_n=o.evals_per_n)

    design = design_engine(inp)
    print_report(design)
    sec = sections_from_table(th)
    wall_t = tuple(s[2] for s in sec)
    start = tuple(s[0] for s in sec) + tuple(s[1] for s in sec)
    base_th = _with_channels(design, th, th.n_channels, start, wall_t)
    base = inj.solve_feed_system(design, base_th, inj.InjectorInputs())["thermal"]

    res = optimize_channels(design, th, inj.InjectorInputs(), wall_t, opt, start_mm=start, progress=print)
    rt.print_thermal_report(design, res["th"], res["feed"]["thermal"])
    inj.print_feed_report(design, res["feed"])
    print_optimize_report(res, dict(t_wg=_peak_t_wg(base), dp_bar=base["dp_regen"] / BAR))


if __name__ == "__main__":
    main()
