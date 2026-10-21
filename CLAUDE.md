# regen_code

Regen-cooled liquid rocket engine design tool. Read README.md first: it covers the architecture,
the models and the known limitations.

- Run everything from this folder with `py -3.14` (plain `python` is 3.11, which lacks CoolProp
  and sv-ttk).
- Non-standard conventions (keep them everywhere): O/F means the CORE O/F, and the fuel film flow is a
  % of the CORE mass flow (core ox + core fuel), e.g. 1 kg/s core + 0.15 kg/s film = 15 %. See README.
- All physics lives in shitty_rocket_code.py, regen_thermal.py, film_cooling.py, injector.py,
  wall_stress.py, materials.py and isopropanol.py. rocket_gui.py only
  collects inputs and displays results; keep physics out of it.
- Coolant models share one interface (regen_thermal.Coolant: props, mu, h, t_from_h, t_sat, plus
  h_liquid and saturation for the film model). film_cooling.py must go through it, never PropsSI by
  name, because isopropanol is a built-in model that CoolProp lacks. Fuel `IPA` is this tool's liquid
  card; RocketCEA's `Isopropanol` card uses the gas-phase heat of formation.
- Inputs are dataclasses (EngineInputs, ThermalInputs, InjectorInputs). Defaults there feed the
  CLIs and the GUI, so CLI argparse defaults must read from the dataclass, not hard-code None.
- GUI colours come from the `THEMES` / `PLOT_COLORS` palettes (light and dark). sv-ttk repaints
  plain tk widgets on every theme change, so create tk widgets through
  `self._reg(widget, bg="surface", fg="muted", ...)` with palette role names, never literal colours.
  Use `PC["..."]` for plot colours, `self._card(..., role="ok")` for tile colours and
  `self._set_status(text, role)` for the status bar.
- Pixel sizes (widths, wraplength) go through `self.px(n)` so they scale with Windows DPI.
- Test GUI changes headlessly: create `App()`, call `withdraw()`, call `run()`, pump `update()`
  until `app.result` is set, then save figures with `fig.savefig`. Never take screen captures. Stub
  `rocket_gui._save_prefs` in tests (it backs `_save_theme` too) so they don't overwrite the user's
  preferences. Layout checks (widget sizes, scrollbars, wheel scrolling) need a mapped window: use
  `app.attributes("-alpha", 0)` instead of `withdraw()`. Only the `__main__` entry point calls
  `restore_window()`, so tests never read or save the user's window state.
- Matplotlib figures are made at `dpi=96`; matplotlib applies the display scaling itself. Don't pass
  `96 * ui_scale` (that double-scales and the plot overflows its canvas).
- After changing physics, check that `py -3.14 shitty_rocket_code.py` output is unchanged and that
  the thermal energy balance still closes.
