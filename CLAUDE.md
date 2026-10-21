# regen_code

Regen-cooled liquid rocket engine design tool. Read README.md first: it covers the architecture,
the models and the known limitations.

- Run everything from this folder with `py -3.14` (plain `python` is 3.11, which lacks CoolProp
  and sv-ttk).
- Non-standard conventions (keep them everywhere): O/F means the CORE O/F, and the fuel film flow is a
  % of the CORE mass flow (core ox + core fuel), e.g. 1 kg/s core + 0.15 kg/s film = 15 %. See README.
- All physics lives in shitty_rocket_code.py, regen_thermal.py, film_cooling.py, injector.py,
  wall_stress.py and materials.py. rocket_gui.py only
  collects inputs and displays results; keep physics out of it.
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
  `rocket_gui._save_theme` in tests so they don't overwrite the user's theme preference.
- After changing physics, check that `py -3.14 shitty_rocket_code.py` output is unchanged and that
  the thermal energy balance still closes.
