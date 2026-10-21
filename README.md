# Regen Engine Designer

Preliminary design tool for a regeneratively cooled liquid rocket engine. From chamber pressure,
thrust, propellants, O/F, L* and efficiencies, it:
- predicts delivered performance
- sizes a conical chamber and nozzle
- runs a steady 1-D regen-cooling thermal circuit
- solves the fuel-side feed pressure and sizes plain-orifice injector holes for fuel and oxidizer
- checks the hot wall for stress, safety factors, low-cycle fatigue life and buckling
- can optimise the cooling channels
- estimates fuel film cooling and sizes the film-cooling holes

> **Non-standard conventions — read before entering O/F or film flow**
> - **O/F is the CORE O/F**: core oxidizer flow / core fuel flow. Film-cooling fuel is *not* included.
> - **Film flow is a percentage of the CORE mass flow** (core oxidizer + core fuel), not of the fuel flow and
>   not of the total flow: ṁ_film = (film % / 100)·ṁ_core. Example: a 1 kg/s core with 0.15 kg/s of film
>   is **15 %**.
> - The core (thrust, throat, ṁ_core, O/F) is sized exactly as without film. The film is extra fuel on top;
>   the reports also give the overall O/F, ṁ_ox / (ṁ_fuel,core + ṁ_film), and the total flows.

## Running it

All commands are run from this folder with Python 3.14.

```
py -3.14 rocket_gui.py                 # desktop GUI (recommended)
py -3.14 shitty_rocket_code.py         # engine only: performance + geometry
py -3.14 regen_thermal.py              # engine + thermal circuit (fixed regen inlet pressure)
py -3.14 injector.py                   # engine + thermal + feed budget + orifice sizing
py -3.14 regen_optimize.py             # channel optimizer (N, width, height) + the injector.py output
py -3.14 wall_stress.py                # injector.py output + hot-wall stresses, safety factors,
                                       #   fatigue life and buckling margin
```

Every command-line script accepts every engine flag (`--pc 30 --of 2.3 --thrust 5000 --film 10 ...`;
see `py -3.14 shitty_rocket_code.py --help`). `--of` is the core O/F and `--film` the film flow in % of
the core flow (see the conventions above). Script-specific flags:

| Script | Flags |
|---|---|
| `regen_thermal.py` | `--channels CSV`, `--coolant-table CSV`, `--thermal-csv`, `--thermal-plot-path` |
| `injector.py` | `--channels CSV`, `--coolant-table CSV` |
| `regen_optimize.py` | `--n-min`, `--n-max`, `--n-step`, `--w-min`, `--h-min`, `--h-max`, `--rib-min`, `--dp-max`, `--opt-stations`, `--evals-per-n`, plus `--channels`, `--coolant-table` |
| `wall_stress.py` | `--material` (e.g. `GRCop-42`), `--t-jacket` (jacket outer temperature, K), `--life-factor`, `--sf-yield`, plus `--channels`, `--coolant-table` |

Defaults come from the dataclasses at the top of each file: `EngineInputs`, `ThermalInputs`,
`InjectorInputs`, `OptimizerInputs`, `StressInputs`. Editing those also changes the GUI defaults.

Dependencies: `rocketcea` (NASA CEA), `CoolProp`, `numpy`, `scipy`, `matplotlib`, and optionally
`sv-ttk` (Windows 11 theme for the GUI; it falls back to the native theme without it).

```
py -3.14 -m pip install rocketcea CoolProp numpy scipy matplotlib sv-ttk
```

## Using the GUI

**Inputs (left)**
- **Input tabs:** Engine, Geometry, Cooling, Channels, Injector. Hover over any input for a
  description and typical range. `auto` lets the code choose the value. The example cases on the
  Engine tab fill in complete inputs.
- **Cooling → Liner wall:** pick the liner material. It fills in the room-temperature conductivity
  (the thermal circuit uses one constant k) and sets the material for the structural checks.
- **Cooling → Wall stress:** tick *Check wall stresses* to run the stress, fatigue and buckling
  checks. Set the yield safety-factor target (used to colour its tile), the jacket outer temperature
  (for the longitudinal stress) and the fatigue life factor.
- **Cooling → Film cooling:** film flow as **% of the core flow** (0 = no film), the film efficiency
  (fraction that stays on the wall) and a gas-to-film h multiplier. The Engine tab's O/F is the core O/F.
- **Injector → Orifices:** the *Film* column sets the number of film-cooling holes and their Cd; they are
  sized like the fuel orifices (same manifold) whenever the film is on.
- **Cooling → Anchoring multipliers:** scale Bartz h_g and coolant h_c to match test data.
- **Channels:** enter width, height and wall thickness at the cylinder, throat and exit. Tick
  *Size widths from the rib thicknesses* to enter a rib per section instead; the widths are then
  computed and written back after each run.
- **Channels → Channel optimizer:** tick *Optimize channels when running* to have Run analysis search
  the channel count and the width and height at each section within the limits given, write the result
  into the table and analyse it. Untick it to use the table as entered (manual or rib-sized widths).

**Outputs (right).** Click **Run analysis**, then read the output tabs:
- **Summary:** performance, geometry, thermal and feed tiles, and the engine, channel, thermal, feed
  and optimizer reports with their warnings. With film on, the Isp tiles show the delivered Isp with
  film (core Isp below it), the *Liquid film length* tile gives the film-covered length and the film
  holes, and the reports add the film flows, liquid film length, film effectiveness and film-hole sizes.
  The plots are unchanged.
- **Structures:** a note that all structural results are fully elastic, then tiles for the governing
  von Mises stress, maximum principal stress, yield safety factor, fatigue life and buckling margin,
  each with its location. Below them is the structural report, which opens with the same elastic
  disclaimer and has these sections:
  1. inputs and load case
  2. definitions of every checked stress
  3. load components at key stations (Δp, σ_b at the rib root and mid-span, σ_t, σ_l, E, α, ν, S_c,
     buckling margin)
  4. σ_x, σ_y, σ_vm, σ_1, σ_3, σ_p, SF_y, strain range and cycles to failure at each of the four
     check points
  5. governing values
  6. low-cycle fatigue
  7. longitudinal buckling
  8. warnings
- **Structural plots** (same elastic note at the top), with three sub-tabs:
  - **Envelope & yield SF:** worst-point σ_vm, σ_p and SF_y along the engine.
  - **Stress components:** each load on its own, and σ_vm at each check point.
  - **Fatigue & buckling:** cycles to failure, allowable firings, and buckling margin.
- **Engine, Temperatures, Heat flux, Coolant, Pressure budget, O/F sweep:** plots. Use the toolbar
  under each to zoom and pan.

**Other controls**
- **Save settings… / Load settings…** write every input (not the results) to a JSON file and read
  it back. Inputs missing from an older file keep their current values; unknown entries are ignored.
- **Export results…** saves the contour and thermal CSVs, the text report (engine, thermal, feed and
  structural) and every plot as PNG.
- **Dark mode** (switch in the header) is the default. The choice is saved to
  `~/.regen_engine_designer.json`.
- The GUI declares Windows DPI awareness and scales its layout and plots, so it stays sharp at
  125 % / 150 % display scaling. Plots are drawn when their tab is first shown after a run.

## How the files fit together

```
shitty_rocket_code.py   EngineInputs -> design_engine() -> EngineDesign (perf, sizing, chamber, nozzle)
        │
regen_thermal.py        ThermalInputs + EngineDesign -> run_thermal() -> station-by-station results
        │               (film_cooling.film_profile() sets the film driving temperature when film > 0)
        │
injector.py             InjectorInputs -> solve_feed_system() -> wraps run_thermal (secant on the
        │                                  regen inlet pressure), returns pressure budget + orifice
        │                                  sizes, and wall stresses when stress=StressInputs(...) is given
        ├── regen_optimize.py   OptimizerInputs -> optimize_channels() -> best N and channel sizes
        │                       (calls run_thermal, then solve_feed_system to verify)
        └── wall_stress.py      StressInputs + thermal result -> wall_stresses() -> stresses, safety
                                factors, fatigue life, buckling margin (material data: materials.py)
rocket_gui.py           Collects inputs, calls design_engine, optimize_channels (if ticked),
                        solve_feed_system and wall_stresses in a worker thread, shows tiles,
                        text reports and plots. No physics.
```

| File | Key functions |
|---|---|
| `shitty_rocket_code.py` | `design_engine`, `evaluate_performance`, `sweep_of`, `chamber_contour`, `nozzle_contour`, `print_report`, `parse_args` |
| `regen_thermal.py` | `run_thermal`, `solve_station`, `bartz_h`, `make_coolant`, `Coolant`, `TableCoolant`, `channel_table_from_sections`, `widths_for_ribs`, `print_thermal_report` |
| `film_cooling.py` | `film_profile`, `hatch_papell_eta`, `velocity_ratio_factor`, `pseudo_critical_temperature` |
| `injector.py` | `solve_feed_system`, `size_orifices`, `ox_density`, `print_feed_report` |
| `regen_optimize.py` | `optimize_channels`, `OptimizerInputs`, `print_optimize_report` |
| `wall_stress.py` | `wall_stresses`, `StressInputs`, `print_stress_report`, `POINTS` (check points and signs), `DEFINITIONS` (report text for each stress) |
| `materials.py` | `MATERIALS`, `get`, `prop` (E, α, ν, yield, ultimate, RA versus temperature; room-temperature k), `cycles_to_failure`, `LCF`, `StrainLife` |
| `isopropanol.py` | `IsopropanolCoolant` (IPA liquid with the `Coolant` interface; CoolProp has no IPA), property functions `density`, `viscosity`, `conductivity`, `cp_liquid`, `enthalpy`, `p_sat`, `t_saturation`, `h_fg` (Appendix B) |
| `rocket_gui.py` | `App` (Tkinter). `_collect` builds the dataclasses; `_work` runs the analysis; `_fill_stats` and `_fill_structures` fill the Summary and Structures tabs |

Coordinates: x = 0 at the throat, negative toward the injector. SI units internally; inputs and
reports use bar, mm, MPa and K.

## Models

**Performance (`shitty_rocket_code.py`).** RocketCEA gives ideal c*, Cf, Tc, γ, MW and Pe in
shifting-equilibrium, frozen, or frozen-at-throat modes. The expansion ratio comes from the design
exit pressure (default: Pe = Pa) or is given directly. Delivered c* = η_c*·c*. Delivered Cf =
η_div·η_nozzle·Cf_vac − ε·Pa/Pc, with η_div = (1 + cos α)/2. Throat area comes from
F/(Pc·Cf), mass flow from Pc·At/c*.

**Geometry.** Conical only, with sharp corners (no fillets or throat radii). Chamber volume =
L*·At, including the convergent frustum. The contraction ratio is an input or comes from Humble's
correlation, Ac/At = 8.0·Dt[cm]^-0.6 + 1.25.

**Thermal circuit (`regen_thermal.py`).** At each axial station:
T_aw → 1/h_g (Bartz) → T_wg → t_w/k_w → T_wc → 1/h_c,eff → T_b.
- Gas side: Bartz with the σ correction, frozen chamber transport properties from CEA, local Mach
  from the area ratio, recovery factor Pr^(1/3). T0 = Tc·η_c*² by default.
- Wall: planar conduction with one constant conductivity k_wall (the material's room-temperature
  value when picked in the GUI).
- Coolant side: Sieder–Tate (default) or Dittus–Boelter in rectangular channels. The ribs are
  treated as straight fins: h_c,eff = h_c(w + 2η_f·h)/(w + L).
- Each station is iterated on T_wg and T_wc. The coolant enthalpy is marched in the flow direction
  (counterflow by default). Pressure drop is Darcy friction with the Haaland friction factor.
- Coolant = the fuel at the full fuel flow. Properties come from CoolProp, auto-mapped from the CEA
  fuel name (RP-1 → n-Dodecane surrogate), from the built-in isopropanol model for IPA (CoolProp has
  none; `isopropanol.py`, Appendix B), or from a CSV table (T_K, rho, mu, cp, k).
- Channel geometry is given at the cylinder (held constant), the throat and the exit, and
  interpolated linearly. Optionally, widths are auto-sized from a rib thickness given at each of
  those sections: w = 2π(r + t_w)/N − rib. The rib then varies linearly between sections.
- Anchoring multipliers `hg_scale` (× Bartz h_g) and `hc_scale` (× coolant h_c, applied before the
  fin correction) default to 1.0. Use them to tune the model to test data or CFD.

**Film cooling (`film_cooling.py`).** A non-entrained, non-reacting fuel film injected at the injector
face, taken from the fuel manifold after the regen jacket (so the jacket carries core fuel + film fuel,
and the film enters at the regen outlet temperature). Two regions, marched along the wall:

1. **Liquid film** after Grisson (*Liquid Film Cooling in Rocket Engines*, AEDC-TR-91-1, 1991), as
   reviewed by Giambelli et al. (EUCASS 2023-648). Assumptions: non-reacting coolant, adiabatic
   film–wall interface, uniform core gas, and **no entrainment mass loss**; only convection from the gas
   (no radiation). The liquid is heated, marched in enthalpy with the coolant model's properties at Pc
   (CoolProp, or the built-in isopropanol model; a CSV property table cannot be used, as it has no
   saturation data), to its
   saturation temperature and then evaporates:
   - heating: d(ṁ_L h_L) = h_g (T_aw − T_L) dA
   - evaporation: dṁ_L = −h_b (T_aw − T_sat) dA / h_fg
   - h_g is the local Bartz coefficient. Evaporation blows mass into the boundary layer, which reduces
     it: h_b = h_g·ln(1+B)/B, with B = c_p,g (T_aw − T_sat)/h_fg·(MW_g/MW_c)^0.6 (Giambelli Eqs. 22, 31).
   - **The liquid film length L_c is where ṁ_L reaches zero.** It is reported from the injector face.
2. **Gaseous film** after Hatch & Papell (NASA TN D-130, 1959, Eq. 12). Assumptions: the film is a
   discrete layer with "no mixing", "in the absence of chemical reactions". The coolant is the evaporated
   fuel at T_sat from L_c on:
   - ln η = −[X − 0.04]·(S·V_g/α_c)^0.125·f(V_g/V_c), with X = ∫h_g dA/(ṁ_c c_p,c) and
     η = (T_aw − T_eff)/(T_aw − T_sat).
   - f = 1 + 0.4·tan⁻¹(V_g/V_c − 1) for V_g/V_c ≥ 1, and (V_c/V_g)^(1.5(V_c/V_g − 1)) for V_g/V_c ≤ 1.
   - The vapour layer is taken to move at the gas velocity (f = 1, Hatch & Papell's optimum), so
     S·V_g/α_c = ṁ_c c_p,c/(2πr k_c) from continuity.
   - Lucas & Golladay (NASA TN D-1988, 1963) found this correlation directly applicable in a small
     JP-4/gaseous-oxygen rocket for a non-reactive coolant with convective heating. They needed a somewhat
     higher h in the convergent section and throat; *film_h_scale* anchors that.
3. **Coupling to the regen circuit**, after Stechman, Oberstone & Howell (*J. Spacecraft Rockets* 6(2),
   1969), whose small-engine model combines the same two regions: q'' = h_g (T_eff − T_wg). The circuit
   keeps Bartz h_g and replaces the gas recovery temperature T_aw with the film driving temperature: T_L
   over the liquid film, then T_eff. Upstream of the film nothing changes. Reports still show the gas T_aw.
- *film_efficiency* (default 1) is the fraction of the film that stays on the wall: an anchoring factor,
  since entrainment is not modelled.
- If Pc is above the fuel's critical pressure (for example the n-dodecane surrogate for RP-1,
  p_crit ≈ 18 bar), the film has no latent heat. The liquid region then ends at the pseudo-critical
  temperature, with a warning that this is an extrapolation.
- Performance: the film adds mass flow but no thrust, so delivered Isp with film = F/(ṁ_core + ṁ_film).
  The throat is sized on the core flow; film blockage of the throat is neglected (a warning says so).

**Feed system and injector (`injector.py`).** Manifold pressure = Pc·(1 + stiffness), with separate
fuel and oxidizer stiffness. The regen inlet pressure is solved by a secant iteration so the fuel
reaches its manifold at that pressure. With film, the same iteration also feeds the regen outlet
temperature back as the film injection temperature until it changes by less than 0.5 K. Orifices are
plain incompressible holes: A = ṁ/(Cd·√(2ρΔP)), split over N holes. Fuel density is taken at the regen
outlet state; oxidizer density comes from CoolProp (LOX, GOX, N2O, F2) or is entered. The fuel orifices
carry the core fuel. Film holes carry the film flow, with their own hole count and Cd, and the fuel
manifold's ΔP and density.

**Channel optimizer (`regen_optimize.py`).** Minimises the peak hot-wall temperature over the
channel count N and the width and height at the cylinder, throat and exit, with wall thickness held
as entered. Limits: minimum width, minimum and maximum height, minimum rib, maximum regen ΔP. The rib
limit is an exact bound on width (the land is linear between sections). For each N, COBYQA (scipy)
runs on the six sizes with ΔP as a constraint, using `run_thermal` at a reduced station count and an
inlet pressure of Pc(1 + s_f) + ΔP_max. N is scanned coarsely, then refined in steps of 1 around the
best. The optimum is re-checked at the full station count with `solve_feed_system`.

**Wall stresses (`wall_stress.py`).** Stresses in the hot wall (liner) at every cooled station, from
the converged thermal solution. Steady firing only. Positive = tension; x is tangential (across the
channel span), y longitudinal, z radial.

> **All structural results assume fully elastic behaviour (σ = Eε), even past yield.** Stresses above
> the yield strength are elastic-equivalent values, not real stresses: the actual stress stays near
> yield and the excess becomes plastic strain. A yield safety factor below 1 means the wall yields, not
> how much. Fatigue life uses the elastic strain range. That is close to the true total strain for
> these thermally driven loads, but it misses plastic strain concentration, ratcheting and creep. The
> buckling margin uses E in place of the tangent modulus, so it is optimistic wherever the wall has
> yielded. The same disclaimer heads the structural report and the Structures and Structural plots
> tabs.

| Stress | Equation | Physical model | Source |
|---|---|---|---|
| Pressure bending σ_b (rib root), σ_b,mid | σ_b = (Δp/2)(w/t_w)², σ_b,mid = σ_b/2 | Wall between two ribs = beam built in at both ends, span w, depth t_w, uniform load Δp. Root moment Δp·w²/12, mid-span Δp·w²/24, σ = 6M/t² | Built-in beam, uniform load (Roark's Formulas for Stress and Strain, Table 8.1) |
| Through-wall thermal σ_t | Eα·q·t_w / (2(1−ν)k) = Eα(T_wg − T_wc) / (2(1−ν)) | Linear temperature drop across a plate restrained against bending; compressive on the hot face, tensile on the cold face | Huzel & Huang, *Design of Liquid Propellant Rocket Engines*, NASA SP-125 (2nd ed. 1971), Eqs. (4-27), (4-31), pp. 107–109 |
| Longitudinal thermal σ_l | Eα·(T_wg − T_jacket) | Gas-side wall's axial growth restrained by the jacket; compressive in the hot wall. T_jacket is the jacket outer temperature input, default ambient (293.15 K) | SP-125 Eq. (4-28), p. 108: ΔT between the gas-side wall (zone I) and the back side and shell (zone II) |
| Tangential σ_x, longitudinal σ_y | σ_x = k_b·σ_b + k_t·σ_t, σ_y = −σ_l, σ_z = 0 | Signed combination at each check point (table below) | — |
| von Mises σ_vm | √(σ_x² − σ_xσ_y + σ_y²) | Distortion-energy criterion with σ_z = 0 at the free faces (thesis Eqs. 9–10 with signs) | von Mises (1913) |
| Principal σ_1, σ_3, σ_p | σ_1 = max(σ_x, σ_y, 0), σ_3 = min(σ_x, σ_y, 0) | No shear on the faces, so σ_x, σ_y and 0 are principal. σ_p = maximum principal stress = the one of largest magnitude (σ_1 or σ_3), with its sign | — |

| Check point | k_b | k_t | Strength taken at |
|---|---|---|---|
| Hot face, rib root | −1 | −1 | T_wg |
| Hot face, mid-span | +½ | −1 | T_wg |
| Cold face, rib root | +1 | +1 | T_wc |
| Cold face, mid-span | −½ | +1 | T_wc |

- Δp = coolant pressure − local gas static pressure, p_gas = Pc(1 + (γ−1)/2·M²)^(−γ/(γ−1)), so the
  pressure load is largest in the nozzle where the gas has expanded.
- In the planar thermal circuit q·t_w/k = T_wg − T_wc exactly (k = the constant room-temperature
  k_wall), so σ_t is computed from the two wall temperatures and stays consistent with the thermal
  solution.
- E, α and ν are taken at the mean wall temperature; strengths at the temperature of the face checked.
  For σ_l, α is the mean CTE from room temperature taken at T_wg, so Eα(T_wg − T_jacket) is the
  restrained thermal strain from an ambient jacket up to the hot-wall temperature.
- **Safety factor.** SF_y = S_yield(T_face) / σ_vm (yield checked against von Mises). Each station
  keeps the worst of its four points. No ultimate-strength factor is reported: in a regen liner the
  elastic stresses exceed yield, where they are not real stresses, so a ratio to the ultimate strength
  would not mean anything. The report and plots give each governing value with its
  location, for example *cold face at rib root, throat, x = 0.0 mm*. For σ_p they also give its
  direction and sense, for example *longitudinal, compressive*.
- **Thesis formulation.** The report also prints the magnitude sum √(½[(|σ_b|+σ_t−σ_l)² + σ_l² +
  (|σ_b|+σ_t)²]) for reference. With signs, the cold face (tangential tension against longitudinal
  compression) usually governs and gives a higher von Mises stress than the magnitude sum.
- Wall material comes from `materials.py`. If `k_wall` differs from the material's room-temperature
  k, a warning says so; the stresses still use the temperatures that were computed.

**Low-cycle fatigue (`wall_stress.py`, `materials.cycles_to_failure`).** One firing is one cycle,
from rest to steady firing and back. At each of the four check points the total strain range is the
elastic equivalent strain Δε = σ_vm / E, and the station life is the shortest of the four. The life
N_f comes from:
- **Measured strain–life curves.**
  - CuCrZr: Coffin–Manson fit Δε/2 = 3.5286·N^−0.262 % (SAA, room temperature to 500 °C; Zhang et
    al. 2019, Eq. 7).
  - GRCop-42: Δε = 16.08·N^−0.3323 %, digitised from the wrought-GRCop-42 line in Gradl et al.
    AIAA 2019-4228, Fig. 2 (room temperature to 600 °C pooled, 561–19 400 cycles). The figure's axis
    reads mm/mm, but the values are percent.
  - AlSi10Mg: modified Basquin–Coffin–Manson fit, ε_a = Δε/2 = (σ′_f/E)(2N)^b + ε′_f(2N)^c + ε_fat, with
    σ′_f = 300 MPa, b = −0.30, ε′_f = 0.10, c = −0.4, E = 70 GPa and the endurance strain ε_fat =
    1.06·10⁻³ (Concli et al., *Int. J. Fatigue* 177 (2023) 107931, Eq. 23, Tables 10–11). As built,
    as-built surface, room temperature only. Below ε_fat the life is unlimited. N is found by bisection.
    See Appendix A.
- **Manson's method of universal slopes** (Manson, *Exp. Mech.* 5 (1965) 193) where only tensile data
  exist: Δε = 3.5(S_u/E)·N^−0.12 + D^0.6·N^−0.6, with ductility D = ln(1/(1−RA)) and S_u, E and RA at
  the face temperature. Inconel 718 uses Special Metals Table 20 RA; 6061-T6 uses Summers et al.
  Fig. 11 RA. This is an estimate and is labelled as one.
- **No fatigue data:** Copper C101 and SS 316L (no reduction-of-area data in the sources); life is
  not computed.

Allowable firings = N_f / life factor (default 4, a scatter factor on mean-curve life). Lives outside
the test data's cycle or temperature range are flagged as extrapolated.

**Longitudinal buckling (`wall_stress.py`).** SP-125 Eq. (4-29) gives the critical stress for
inelastic buckling of the gas-side wall, S_c = 4E_tE_c·t / ((√E_t + √E_c)²·√(3(1−ν²))·r). It is
printed in SP-125 with a typo, (√E_t √E_c)², which makes it dimensionally inconsistent. With no
tangent moduli in the library, E_t = E_c = E, so S_c = E·t/(r·√(3(1−ν²))). That is classical elastic
axial buckling of a cylinder, evaluated with E at the mean wall temperature, wall thickness t_w and
mean radius r + t_w/2. SP-125 keeps σ_l ≤ 0.9·S_c, so the margin is MS_b = 0.9·S_c/σ_l − 1. Wherever
σ_l exceeds yield the tangent modulus collapses and the real critical stress is much lower; the
elastic margin is then an unconservative upper bound, and a warning says so.

**Materials (`materials.py`).** Room-temperature conductivity (the thermal circuit uses one constant
k) and, versus temperature: E, mean CTE α from room temperature, ν, 0.2 % yield, ultimate tensile
strength (used only by the universal-slopes fatigue estimate) and, where available, reduction of area. Lookups are linear in temperature. Outside a
table's range the end value is used, with a warning. A property known at only one temperature is held
constant, also with a warning.

| Material | Condition | k [W/m·K] | Strength data | Fatigue | Sources |
|---|---|---|---|---|---|
| CuCrZr | Solution annealed + aged (ITER grade) | 320 | 20–700 °C | Measured (Coffin–Manson) | E, yield, UTS, fatigue: average fits, Zhang et al., *Fusion Eng. Des.* 144 (2019) 148, Eqs. (1), (3), (5), (7). k, α, ν: Bäker et al., *Mech. Adv. Mater. Mod. Process.* 1:5 (2015), Table 1 |
| GRCop-42 | Wrought; L-PBF reported equivalent | 346 | 20–800 °C | Measured (digitised line) | Yield, UTS, α, k, fatigue: Gradl et al., AIAA 2019-4228, Figs. 1, 2, 5, 6, 7, digitised from the vector figures. E: Mahtabi et al., *JOM* 77 (2025), room temperature, held constant. ν: pure copper (NIST Monograph 177) |
| Copper C101 | Annealed OFHC | 390 | 20–482 °C | None | Yield 34 MPa (5 ksi) and UTS 69 MPa (10 ksi, held constant): NASA Tech Brief 66-10273 (1966), braze-coarsened OFHC. E, α, ν at 295 K: NIST Monograph 177. k: typical handbook value, not re-verified |
| Al 6061-T6 | T6 / T651 | 174 | 24–427 °C | Universal slopes (estimate) | Yield, UTS: Kaufman (2000) data in Summers et al., *Fire Sci. Rev.* 4:3 (2015), Figs. 9–10, read from the figures. RA: Fig. 11. E: Eurocode 9 values in the same paper, Fig. 8. α: EN 1999-1-2 thermal-elongation equation. ν = 0.3: EN 1999-1-1. k: Summers Fig. 1 |
| AlSi10Mg | L-PBF, as built | 120 | 25–400 °C | Measured (Basquin–Coffin–Manson, room temperature) | E, yield, UTS: Cao et al., *J. Mater. Sci. Technol.* 62 (2021) 162, Table 3. k: Martínez-Maradiaga et al., *J. Mater. Eng. Perform.* 29 (2020) 7125, Table 1. Fatigue: Concli et al., *Int. J. Fatigue* 177 (2023) 107931. α: EOS process data sheet (M 290, 30 µm). ν = 0.3: EN 1999-1-1. Details in Appendix A |
| Inconel 718 | Annealed 1800 °F + aged | 11.4 | 27–704 °C | Universal slopes (estimate) | Special Metals SMC-045 (2007): Tables 4 (E, ν), 5 (k, α), 20 (yield, UTS, RA) |
| SS 316L | Annealed, typical 316/316L | 16.2 (at 100 °C) | 20–871 °C | None | AK Steel 316/316L bulletin: Tables 1–2 (yield, UTS), α, k. E 193 GPa at room temperature, held constant. ν = E/2G − 1 |

Values read from raster figures are good to about ±5 MPa (strength) or ±1.5 GPa (E). Values held
constant (GRCop-42 and SS 316L E, all C101 elastic and thermal values) over-predict thermal stress
where E or α would fall with temperature, and under-predict it where α would rise.

## Known limitations

- Sharp throat: Bartz's curvature term uses a nominal `throat_curv_ratio`, so the throat heat-flux
  peak is sensitive to it.
- CoolProp's n-Dodecane model stops at 700 K. Wall-side property lookups above that are clamped,
  with a warning. Supply a real RP-1 property table for better accuracy.
- Not modeled: radiation, axial conduction, wall curvature, temperature-dependent wall conductivity,
  boiling and coking, manifold or feed-line losses, and compressible (gas) injector
  flow (a warning is shown).
- Film cooling is a preliminary estimate. It has no entrainment, no radiation, no mixing or reaction of
  the vapour with the core, and no film thrust. The vapour layer is assumed to move at the gas velocity.
  Liquid properties and saturation are taken at Pc. Supercritical films are extrapolated. The film's own
  heat budget assumes an adiabatic film–wall interface (Grisson), so heat exchanged between the film and
  the regen-cooled wall is not returned to the film. The transpiration (blowing) correction roughly
  doubles the liquid film length in the test case, so anchor film_efficiency and film_h_scale to data
  where possible.
- Wall stresses are fully elastic (even past yield) and closed-form, for steady firing only. Not
  included: plasticity,
  ratcheting or creep; stress concentration at the rib fillets; closeout, rib or jacket stresses;
  start-up or shutdown transients. The jacket is treated as rigid at its outer temperature, which
  makes σ_l conservative. Copper liners usually yield in compression on the hot face each firing, so
  SF_y < 1 there calls for a look at the fatigue life, not necessarily a redesign.
- Fatigue life uses the elastic strain range (no plastic redistribution, mean-stress or hold-time
  effect, no ratcheting or wall thinning), all known to shorten copper liner life. The universal-slopes
  lives are estimates.
- Buckling uses the elastic modulus, so the margin is optimistic wherever σ_l exceeds yield.
- Material data are typical or average values from the sources above, not design minimums
  (A- or B-basis). Additively made properties depend strongly on process, orientation and heat
  treatment.
- AlSi10Mg fatigue data are at room temperature only, so lives at hot points are flagged as extrapolated.
  Its solidus is about 830 K (≈ 557 °C), and its strength is nearly gone by 400 °C, so the hot-wall
  limit for an AlSi10Mg liner should be set far below the 800 K default.
- Isopropanol coolant (Appendix B): saturated-liquid properties with no pressure correction. Liquid c_p
  above 473 K is extrapolated, and the model is capped at 0.98 T_c (498 K). There is no
  supercritical/pseudo-critical model, and boiling is flagged but not modelled. IPA's critical point
  (508 K, 47.6 bar) is close to typical regen outlet conditions, so check the boiling warnings.

## Checks that have passed

- Thermal energy balance: Σq''·dA matches ṁ_c·Δh within 0.03 %.
- Orifice area and auto-sized channel widths agree with hand calculations. Per-section ribs come
  out exactly at the cylinder, throat and exit and vary linearly between them.
- Regen outlet pressure matches the injector target to within 1e-3 bar.
- Adding the structural and film-cooling modules left the output of `shitty_rocket_code.py`,
  `regen_thermal.py` and `injector.py` byte-identical with the film off.
- Film cooling:
  - Hatch & Papell's worked example (TN D-130) is reproduced: η = 0.7212 against 0.722, and coolant
    flow 16.87 against 16.84 lb/s.
  - The liquid film's absorbed heat equals ṁ_film·(Δh_sensible + h_fg) exactly.
  - With blowing switched off, the marched liquid length matches Stechman's closed form: 66.0 against
    66.5 mm for LOX/ethanol, 2 kN, 10 % film.
  - Film-hole area matches a hand calculation: 5.93 mm².
- Wall stresses at the throat (CuCrZr, default case) agree with an independent hand calculation from
  the source equations at all four check points: Δp, σ_b, σ_t, σ_l, von Mises, σ_p and the yield
  safety factor within 0.01 %. The thesis Eq. 10 value matches too.
- Fatigue life and buckling stress agree with hand calculations. CuCrZr at Δε = 0.874 % gives
  N ≈ 2 900 cycles from the Coffin–Manson fit. S_c in the chamber implies E = 70.4 GPa, the UKAEA fit at
  the 480 °C mean wall temperature.
- AlSi10Mg fatigue: the implemented fit gives about 990 cycles at a strain amplitude of 6.3·10⁻³. Concli
  et al. compare this against other AlSi10Mg data, which show about 1 000 cycles at 4–6·10⁻³.
- Isopropanol: see the checks in Appendix B (handbook spot values; the analytic enthalpy matches
  numerical integration of c_p; the thermal energy balance closes exactly with IPA as the coolant, film
  off and on). With the IPA, AlSi10Mg and film-refactor changes, the default outputs of
  `shitty_rocket_code.py`, `regen_thermal.py` and `injector.py` are unchanged.

## Appendix A: AlSi10Mg material data review

Scope: L-PBF AlSi10Mg in the **as-built** condition, the state a printed liner is in unless it is heat
treated. All values are typical or mean values, not design allowables.

| Property | Before | Now | Source and notes |
|---|---|---|---|
| Thermal conductivity k (room temperature, used by the thermal circuit) | 100 W/m·K (EOS data sheet, vertical) | **120 W/m·K** | Martínez-Maradiaga, Mishin & Engelbrecht, "Thermal properties of selectively laser melted AlSi10Mg products with different densities", *J. Mater. Eng. Perform.* 29 (2020) 7125–7130, Table 1: laser flash measurements, as built, 25 °C. 99.5 % dense: 120 (⊥ build) and 122 (∥ build) W/m·K. 99 % dense: 115 and 115. The paper notes these match die-cast alloys of similar composition. EOS quotes 100–110 W/m·K as manufactured. |
| Strength at temperature (yield, UTS, E) | Cao et al. 2021, Table 3 | **Unchanged** | Cao et al., *J. Mater. Sci. Technol.* 62 (2021) 162–172, Table 3: as built, tested at 25–400 °C after a 30 min hold. Yield 322 / 298 / 236 / 143 / 25 MPa and UTS 460 / 382 / 266 / 150 / 30 MPa at 25 / 100 / 200 / 300 / 400 °C. E 76 → 18 GPa. The trend agrees with Uzan et al., *Addit. Manuf.* 24 (2018) 257 (stress relieved, 25–400 °C): yield falls sharply above 200 °C, and UTS falls steadily. |
| Fatigue | None (life not computed) | **Measured strain–life curve** | Concli, Gerosa, Panzeri & Fraccaroli, "High and low cycle fatigue properties of selective laser melted AISI 316L and AlSi10Mg", *Int. J. Fatigue* 177 (2023) 107931 (open access). Fully reversed LCF (ASTM E606) and HCF, as built with as-built surfaces, room temperature. Modified Basquin–Coffin–Manson (Eq. 23, Table 11): σ′_f = 300 MPa, b = −0.30, ε′_f = 0.10, c = −0.4, ε_fat = 1.06·10⁻³ (the 123 MPa fatigue limit at R = −1 over E = 70 GPa; Tables 6 and 10). Test lives 13 to 8.5·10⁵ cycles. |
| α, ν | EOS; EN 1999-1-1 | Unchanged | — |

Notes for use:
- **Conductivity rises once the wall anneals.** Above about 250 °C the as-built cellular Si network
  breaks up (Martínez-Maradiaga et al., Sec. 4). After one heating to 400 °C, k at 25 °C is 140–152
  W/m·K, and at 400 °C it is 144–189 W/m·K. A liner that has fired is therefore more conductive than
  the 120 used here (conservative for wall temperature), but also weaker: strength drops with the
  same coarsening.
- **The fatigue curve is conservative for a machined or polished channel**: it uses as-built rough
  surfaces and sub-surface porosity. It has no temperature dependence. Hot dwell-fatigue tests on
  L-PBF AlSi10Mg report shorter lives at 250 °C than at 100 °C (cyclic softening), so lives at hot
  points are unconservative as well as extrapolated, and are flagged.
- **Not used:** Khan, Yadav, Gradl, Shao & Shamsaei, *Mater. Des.* 265 (2026) 115902, which has
  strain-controlled fatigue of L-PBF AlSi10Mg from −195 to 200 °C (NASA MSFC co-author) and would be the
  best source for a temperature-dependent curve. Its full text could not be retrieved, and the search
  snippets disagreed on its fit constants, so none were used. Worth adding if you have access.
- The solidus is about 557 °C (830 K). Set `t_wall_limit_K` well below that for an aluminium liner.

## Appendix B: Isopropyl alcohol (IPA) as propellant and coolant

**As a reactant (`shitty_rocket_code.py`).** RocketCEA ships an "Isopropanol" card, but its enthalpy
is −65 133 cal/mol = −272.5 kJ/mol, the **gas-phase** heat of formation. That overstates the energy of
a liquid propellant by the heat of vaporisation. The tool therefore registers its own fuel card,
**`IPA`** (`add_new_fuel`): C₃H₈O, h = −318.1 kJ/mol (liquid, CRC Handbook) = −76 028.7 cal/mol at
298.15 K, ρ = 0.781 g/cm³. No surrogate is needed; CEA has the species. Effect at Pc = 20 bar, ε = 4,
O/F 1.4–2.0:

| O/F | Isp, stock card [s] | Isp, `IPA` card [s] | T_c, stock [K] | T_c, `IPA` [K] |
|---|---|---|---|---|
| 1.4 | 283.1 | 277.2 | 3027 | 2935 |
| 1.6 | 289.2 | 284.8 | 3234 | 3175 |
| 1.8 | 290.5 | 287.2 | 3333 | 3292 |
| 2.0 | 288.4 | 285.8 | 3372 | 3339 |

Select `IPA` (GUI fuel list, or `--fuel IPA`). `Isopropanol` still works but uses the stock card.

**As a coolant (`isopropanol.py`).** CoolProp has no IPA. `IsopropanolCoolant` implements the same
interface as the CoolProp wrapper (props, h, T(h), T_sat, and saturation data for the film model) from
published correlations. Coolant "auto" with fuel `IPA` or `Isopropanol` selects it, as does the GUI
entry "Isopropanol (built-in model)". Coefficients (CAS 67-63-0) were transcribed from the data tables
in the `chemicals` Python package v1.5.2 (C. Bell et al.), which reproduce the primary sources below.
τ = 1 − T/T_c, with T_c = 508.25 K and p_c = 47.62 bar.

| Property | Correlation | Source | Range used |
|---|---|---|---|
| Liquid density | ρ = ρ_c + Aτ^0.35 + Bτ^(2/3) + Cτ + Dτ^(4/3); ρ_c = 273, A = 865.9165, B = −744.0943, C = 975.9052, D = −381.1278 kg/m³ | VDI Heat Atlas, 2nd ed. (Springer 2010), Sec. D3, PPDS | to T_c |
| Liquid viscosity | μ = E·exp[A z^(1/3) + B z^(4/3)], z = (C−T)/(T−D); A = 3.07172, B = 3.84587, C = 580.414 K, D = −62.691 K, E = 7.42·10⁻⁶ Pa·s | VDI Heat Atlas, PPDS eq. 9 | 185–498 K |
| Liquid conductivity | k = 0.2028 − 2.238·10⁻⁴T − 3.1·10⁻⁸T² + 8.3·10⁻¹¹T³ − 6.9·10⁻¹⁴T⁴ | VDI Heat Atlas, PPDS | 185–498 K |
| Liquid c_p, enthalpy | c_p/R = Σ A_j (T/100 K)^(j−1), isobaric cubic splines 185–270, 270–360, 360–473.2 K; h = analytical integral, zero at 298.15 K | Zábranský, Růžička, Majer & Domalski, *Heat Capacity of Liquids: Critical Review and Recommended Values*, J. Phys. Chem. Ref. Data Monograph 6 (1996) | data to 473.2 K, extrapolated to 498 K with a warning |
| Vapour pressure | ln(p/p_c) = (Aτ + Bτ^1.5 + Cτ^2.5 + Dτ^5)·T_c/T; A = −8.44737, B = 1.17402, C = −6.97876, D = 0.69248 | VDI Heat Atlas, PPDS (Wagner 2.5–5 form) | to T_c |
| Latent heat | h_fg = R T_c (Aτ^(1/3) + Bτ^(2/3) + Cτ + Dτ² + Eτ⁶)/M; A = 14.13231, B = −18.02015, C = 33.713095, D = −20.808299, E = −7.959708 | VDI Heat Atlas, PPDS eq. 12 | to T_c |
| Vapour c_p (film model) | ideal gas, c_p/R = 3.334 + 0.018853T + 3.644·10⁻⁵T² − 6.115·10⁻⁸T³ + 2.543·10⁻¹¹T⁴ | Poling, Prausnitz & O'Connell, *Properties of Gases and Liquids*, 5th ed. (2001), App. A | 50–1000 K |
| Vapour conductivity (film model) | k = 0.001447 + 1.1·10⁻⁵T + 1.2455·10⁻⁷T² | VDI Heat Atlas, PPDS (low pressure) | — |

Checks (model value vs handbook value):
- Density at 298 K: 781.8 vs 781 kg/m³ (CRC).
- Viscosity at 298 K: 2.01 vs 2.04 mPa·s.
- c_p at 298 K: 2.58 vs 2.60 kJ/kg·K (CRC 156.5 J/mol·K).
- Conductivity at 298 K: 0.135 vs 0.135 W/m·K.
- p_sat at 355.36 K: 1.013 bar, so T_b = 355.4 K (CRC 355.4 K).
- h_fg at T_b: 677 vs 663 kJ/kg (CRC 39.85 kJ/mol), +2 %.
- The c_p splines are continuous at 270 and 360 K.
- The analytic enthalpy matches trapezoidal integration of c_p to 0.1 J/kg; T(h(T)) returns T.
- LOX/IPA, 4 kN, 20 bar: thermal energy balance Σq''dA = ṁ·Δh closes exactly (490.763 kW film off;
  135.521 kW with 10 % film).

Limitations:
- Pressure effects on the compressed liquid are neglected: saturated-liquid values are used at all
  pressures. This effect has not been quantified here.
- Above 473 K, c_p is an extrapolation. There is no pseudo-critical model above p_c: the film's liquid
  region then ends at 498 K, with a warning.
- IPA saturates at 460 K at 20 bar and 481 K at 30 bar, so small LOX/IPA engines easily reach bulk
  boiling. The existing boiling warnings apply.

## Appendix C: thermal model compared with RPA (no film cooling)

Full report: https://claude.ai/artifact/KHkSS3Ts4YqnMZDjRc9Ha9. RPA (Ponomarenko, *RPA: Tool for Liquid
Propellant Rocket Engine Analysis*, 2012, and the RPA thermal-analysis documentation) uses the same
series thermal circuit and fin form as this tool. The differences that matter, with film cooling
excluded:

- **Gas temperature:** RPA uses the ideal T_c; this tool's default is η_c*²·T_c.
- **Gas-side correlation:** RPA offers Bartz, Ievlev or their average, plus optional H₂O/CO₂ radiation.
- **Coolant-side correlation:** RPA uses fluid-specific Nusselt correlations. For kerosene,
  0.021Re^0.8Pr^0.4(0.64 + 0.36T_c/T_wc); for other fluids, Dittus–Boelter. This tool defaults to
  Sieder–Tate.

Changing these one at a time in this tool (CuCrZr, k = 320 W/m·K) gives the following peak T_wg:

| Case | Default | Ideal T_c | Dittus–Boelter | Both (RPA-like, Bartz mode) |
|---|---|---|---|---|
| LOX/RP-1 1.8 kN, 20 bar | 751 K | 801 K | 874 K | 943 K |
| LOX/ethanol 2 kN, 20 bar | 853 K | 929 K | 1046 K | 1139 K |

For the same engine input, RPA in Bartz mode is therefore expected to predict the **hotter** wall,
mainly through the lower coolant-side coefficient and the higher gas temperature. Radiation, if
enabled in RPA, adds to that. RPA's Ievlev and average modes were not reproduced here, so their
direction relative to this tool is not established.
