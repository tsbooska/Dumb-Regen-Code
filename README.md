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
| `materials.py` | `MATERIALS`, `get`, `prop` (E, α, ν, yield, ultimate, RA versus temperature; room-temperature k), `cycles_to_failure`, `LCF` |
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
  fuel name (RP-1 → n-Dodecane surrogate), or from a CSV table (T_K, rho, mu, cp, k).
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
   (no radiation). The liquid is heated, marched in enthalpy with CoolProp properties at Pc, to its
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
- **Manson's method of universal slopes** (Manson, *Exp. Mech.* 5 (1965) 193) where only tensile data
  exist: Δε = 3.5(S_u/E)·N^−0.12 + D^0.6·N^−0.6, with ductility D = ln(1/(1−RA)) and S_u, E and RA at
  the face temperature. Inconel 718 uses Special Metals Table 20 RA; 6061-T6 uses Summers et al.
  Fig. 11 RA. This is an estimate and is labelled as one.
- **No fatigue data:** Copper C101, SS 316L and AlSi10Mg (no reduction-of-area data in the sources);
  life is not computed.

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
| AlSi10Mg | L-PBF, as built | 100 | 25–400 °C | None | E, yield, UTS: Cao et al., *J. Mater. Sci. Technol.* 62 (2021) 162, Table 3. k, α: EOS process data sheet (M 290, 30 µm). ν = 0.3: EN 1999-1-1 |
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
  boiling and coking, film cooling, manifold or feed-line losses, and compressible (gas) injector
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
