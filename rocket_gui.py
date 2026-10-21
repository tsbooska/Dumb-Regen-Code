"""
Regen Engine Designer - desktop front end.

Front end for shitty_rocket_code.py (performance + conical geometry), regen_thermal.py
(steady regen-cooling thermal circuit) and injector.py (feed pressure budget + orifice sizing).
All physics lives in those files; this one collects inputs, runs them and shows the results.

    py -3.14 rocket_gui.py

Defaults come from EngineInputs / ThermalInputs / InjectorInputs, so editing those dataclasses
changes the GUI defaults.
"""

from __future__ import annotations

import contextlib
import ctypes
import io
import json
import logging
import os
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import injector as inj
import materials as mat
import regen_optimize as ro
import regen_thermal as rt
import shitty_rocket_code as eng
import wall_stress as ws

try:  # Windows 11 style theme (pip install sv-ttk); falls back to the native theme without it
    import sv_ttk
except ImportError:
    sv_ttk = None

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# --------------------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------------------
UI = "Segoe UI Variable Text" if os.name == "nt" else "Segoe UI"
UI_DISPLAY = "Segoe UI Variable Display" if os.name == "nt" else "Segoe UI"
F = (UI, 10)
FB = (UI, 10, "bold")
MONO = ("Cascadia Mono", 10)

# Windows 11 light / dark palettes. Every value must be unique within a palette: switching themes
# recolours existing widgets by mapping each old colour to the new one.
THEMES = {
    "light": dict(bg="#f3f3f3", surface="#ffffff", border="#e5e5e5", fg="#1b1b1b", muted="#616161",
                  accent="#005fb8", ok="#0f7b0f", warn="#9d5d00", bad="#c42b1c"),
    "dark": dict(bg="#202020", surface="#2b2b2b", border="#3a3a3a", fg="#f0f0f0", muted="#a6a6a6",
                 accent="#60cdff", ok="#6ccb5f", warn="#fce100", bad="#ff99a4"),
}
# Plot series keep their hue in both themes but are lifted for contrast on the dark background
PLOT_COLORS = {
    "light": dict(grid="#e5e7eb", edge="#9ca3af", hot="#c0392b", warm="#e67e22", cool="#2471a3", aw="#9ca3af",
                  violet="#6c3483", green="#1e8449", sky="#2e86c1", gas="#fde2c8", liner="#a0522d",
                  channel="#7fb3e0", bar_pc="#d1d5db", bar_inj="#f0a35e", bar_regen="#5b8fd6"),
    "dark": dict(grid="#383838", edge="#6b6b6b", hot="#ff6b5b", warm="#ffad5c", cool="#5fb4f5", aw="#8a8a8a",
                 violet="#c39bf0", green="#5fd38d", sky="#7fd3f5", gas="#4a3324", liner="#c27a4a",
                 channel="#2f6fa8", bar_pc="#4a4a4a", bar_inj="#c98642", bar_regen="#4f86cc"),
}
PREFS_PATH = os.path.join(os.path.expanduser("~"), ".regen_engine_designer.json")

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.linewidth": 0.8,
    "legend.frameon": False, "legend.fontsize": 8,
})


def _use_palette(name: str) -> None:
    """Point the colour globals (read at call time by widgets and plots) and matplotlib at a theme."""
    global THEME, BG, SURFACE, BORDER, FG, MUTED, ACCENT, OK, WARN, BAD, TIP_BG, PC
    t = THEMES[name]
    THEME = name
    BG, SURFACE, BORDER, FG, MUTED = t["bg"], t["surface"], t["border"], t["fg"], t["muted"]
    ACCENT, OK, WARN, BAD = t["accent"], t["ok"], t["warn"], t["bad"]
    TIP_BG = SURFACE
    PC = PLOT_COLORS[name]
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "grid.color": PC["grid"],
        "axes.edgecolor": PC["edge"], "text.color": FG, "axes.labelcolor": FG, "axes.titlecolor": FG,
        "xtick.color": MUTED, "ytick.color": MUTED, "legend.labelcolor": FG,
    })


def _load_theme() -> str:
    try:
        with open(PREFS_PATH, encoding="utf-8") as f:
            name = json.load(f).get("theme", "dark")
    except (OSError, ValueError, AttributeError):
        return "dark"
    return name if name in THEMES else "dark"


def _save_theme(name: str) -> None:
    try:
        with open(PREFS_PATH, "w", encoding="utf-8") as f:
            json.dump({"theme": name}, f)
    except OSError:
        pass


_use_palette(_load_theme())

OXIDIZERS = ["LOX", "GOX", "N2O", "N2O4", "H2O2", "IRFNA", "MON25", "F2"]
FUELS = ["RP-1", "Ethanol", "IPA", "Methanol", "CH4", "LCH4_NASA", "LH2", "GH2", "MMH", "UDMH", "N2H4", "A50",
         "JetA", "Propane"]
COOLANTS = {"Auto (from fuel)": "auto", "RP-1 (n-Dodecane surrogate)": "n-Dodecane", "Ethanol": "Ethanol",
            "Isopropanol (built-in model)": "Isopropanol",
            "Methanol": "Methanol", "Methane": "Methane", "Propane": "Propane", "Hydrogen": "Hydrogen",
            "Ammonia": "Ammonia", "Nitrous oxide": "NitrousOxide", "Oxygen": "Oxygen"}
MATERIALS = {m.label: m.k for m in mat.MATERIALS.values()}   # label -> room-temperature k [W/m-K]

SETTINGS_VERSION = 1  # format of the Save / Load settings JSON

PRESETS = {
    "Defaults": {},
    "LOX / Ethanol, 1 kN": dict(
        ox="LOX", fuel="Ethanol", pc="15", thrust="1000", of="1.5", lstar="1.2", cr="auto",
        coolant="Auto (from fuel)", t_in="290", n_ch="24", n_fuel="12", n_ox="12",
        cyl_w="1.5", cyl_h="2", cyl_t="1", thr_w="1", thr_h="1.2", thr_t="0.8",
        exit_w="1.5", exit_h="2", exit_t="1"),
    "LOX / Methane, 25 kN": dict(
        ox="LOX", fuel="CH4", pc="40", thrust="25000", of="3.2", lstar="1.0", cr="auto",
        coolant="Auto (from fuel)", t_in="120", n_ch="60", n_fuel="60", n_ox="60", t_ox="90",
        cyl_w="2", cyl_h="3", cyl_t="1.2", thr_w="1.5", thr_h="2", thr_t="0.9",
        exit_w="2", exit_h="3", exit_t="1.2"),
}

CARD_NAMES = ["Isp (design)", "Isp (vacuum)", "Thrust", "Mass flow", "c*", "Throat diameter",
              "Exit diameter", "Max hot-wall T", "Regen ΔP", "Feed pressure", "Orifice diameter",
              "Liquid film length"]
STRUCT_CARD_NAMES = ["Max von Mises stress", "Max principal stress", "Min yield SF", "Fatigue life",
                     "Buckling margin"]
ELASTIC_NOTE = ("All structural results assume fully elastic behaviour, even past yield: stresses above the "
                "yield strength are elastic-equivalent values, not real stresses.")


def _enable_hidpi() -> None:
    """Declare DPI awareness so Windows renders the window at native resolution instead of
    bitmap-stretching a 96-DPI image (the cause of blurry Tk apps on scaled displays)."""
    if os.name != "nt":
        return
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(1)  # system DPI aware
    except (AttributeError, OSError):
        try:
            ctypes.windll.user32.SetProcessDPIAware()
        except (AttributeError, OSError):
            pass


def _fmt(v) -> str:
    if v is None:
        return "auto"
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float)):
        return f"{v:g}"
    return str(v)


# --------------------------------------------------------------------------------------
# Small widgets
# --------------------------------------------------------------------------------------
class Tooltip:
    def __init__(self, widget, text: str):
        self.widget, self.text, self.tip = widget, text, None
        widget.bind("<Enter>", self.show, add="+")
        widget.bind("<Leave>", self.hide, add="+")

    def show(self, _=None):
        if self.tip:
            return
        x = self.widget.winfo_rootx() + 16
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.geometry(f"+{x}+{y}")
        scale = self.widget.winfo_fpixels("1i") / 96
        tk.Label(self.tip, text=self.text, bg=TIP_BG, fg=FG, font=(UI, 9), bd=0, highlightthickness=1,
                 highlightbackground=BORDER, wraplength=int(340 * scale), justify="left", padx=8, pady=6).pack()

    def hide(self, _=None):
        if self.tip:
            self.tip.destroy()
            self.tip = None


class ScrollFrame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, bg=BG, highlightthickness=0)
        sb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(win, width=e.width))
        self.canvas.configure(yscrollcommand=sb.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")


# --------------------------------------------------------------------------------------
# The app
# --------------------------------------------------------------------------------------
class App(tk.Tk):
    def __init__(self):
        _enable_hidpi()
        super().__init__()
        self.title("Regen Engine Designer")
        self.ui_scale = self.winfo_fpixels("1i") / 96  # 1.0 at 100 % Windows scaling, 1.5 at 150 %
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry(f"{min(self.px(1500), int(sw * 0.95))}x{min(self.px(920), int(sh * 0.88))}")
        self.minsize(min(self.px(1200), int(sw * 0.9)), min(self.px(720), int(sh * 0.8)))
        self.option_add("*TCombobox*Listbox.font", F)
        self.protocol("WM_DELETE_WINDOW", self._close)

        self.vars: dict[str, tk.Variable] = {}
        self.labels: dict[str, str] = {}
        self.defaults: dict[str, object] = {}
        self.scrollers: list[ScrollFrame] = []
        self.figs: dict[str, tuple[Figure, FigureCanvasTkAgg]] = {}
        self.fig_tabs: dict[str, tuple] = {}  # plot key -> (notebook, its tab's widget name, outer tab or None)
        self.texts: list[tk.Text] = []        # report panes (Summary, Structures)
        self.stale: set[str] = set()         # plots changed since last drawn (drawn when their tab shows)
        self.result = None
        self.jobs: queue.Queue = queue.Queue()
        self.busy = False
        self.toolbars: list[NavigationToolbar2Tk] = []
        self._themed: list[tuple[tk.Misc, dict[str, str]]] = []
        self.dark_var = tk.BooleanVar(value=THEME == "dark")

        self._style()
        self._header()
        self.status = self._reg(tk.Label(self, text="Ready. Hover over any input for a description.",
                                         font=(UI, 9), anchor="w", padx=18, pady=6, bd=0), bg="bg", fg="muted")
        self.status.pack(fill="x", side="bottom")

        body = ttk.Frame(self, padding=(14, 4, 14, 0))
        body.pack(fill="both", expand=True)
        self.left = ttk.Frame(body, width=self.px(470))
        self.left.pack(side="left", fill="y")
        self.left.pack_propagate(False)
        right = ttk.Frame(body)
        right.pack(side="left", fill="both", expand=True, padx=(10, 0))

        self._build_inputs(self.left)
        self._build_buttons(self.left)
        self._build_outputs(right)
        self.bind_all("<MouseWheel>", self._wheel)
        self.after_idle(self._paint)  # after sv-ttk's first repaint of tk widgets
        self._titlebar()

    def px(self, n: float) -> int:
        """Convert a layout size given at 96 DPI to real pixels for this display."""
        return int(round(n * self.ui_scale))

    # ---------------------------------------------------------------- style / header
    def _style(self):
        s = ttk.Style(self)
        if sv_ttk is not None:
            sv_ttk.set_theme(THEME)
        else:
            s.theme_use("vista" if "vista" in s.theme_names() else "clam")
        s.configure(".", font=F)
        s.configure("TFrame", background=BG)
        s.configure("TLabel", background=BG, foreground=FG)
        s.configure("Muted.TLabel", background=BG, foreground=MUTED, font=(UI, 9))
        s.configure("Head.TLabel", background=BG, foreground=MUTED, font=(UI, 9, "bold"))
        s.configure("TLabelframe", background=BG)
        s.configure("TLabelframe.Label", background=BG, foreground=FG, font=FB)
        s.configure("TCheckbutton", background=BG)
        s.configure("TRadiobutton", background=BG)
        if sv_ttk is not None:
            # sv-ttk draws the notebook pane with an image element that Tk tiles over the whole pane on
            # every tab switch (100-400 ms at high DPI). A plain client area keeps its tabs, drops that cost.
            s.layout("TNotebook", [("Notebook.client", {"sticky": "nswe"})])
            # The client element draws a frame in these colours; unset, Tk uses light defaults.
            s.configure("TNotebook", bordercolor=BG, lightcolor=BG, darkcolor=BG, padding=0)
        s.configure("TNotebook", background=BG)
        s.configure("TNotebook.Tab", padding=[12, 4])
        s.configure("Run.TButton", font=(UI, 11, "bold"), padding=(10, 8))
        self.configure(bg=BG)

    def _header(self):
        top = self._reg(tk.Frame(self), bg="bg")
        top.pack(fill="x")
        txt = self._reg(tk.Frame(top), bg="bg")
        txt.pack(side="left", padx=18, pady=(12, 6))
        self._reg(tk.Label(txt, text="Regen Engine Designer", font=(UI_DISPLAY, 18, "bold")),
                  bg="bg", fg="fg").pack(anchor="w")
        self._reg(tk.Label(txt, text="Performance  ·  Geometry  ·  Regenerative cooling  ·  Feed system & injector",
                           font=(UI, 10)), bg="bg", fg="muted").pack(anchor="w")
        ttk.Button(top, text="Help", command=self._help).pack(side="right", padx=(6, 18))
        ttk.Checkbutton(top, text="Dark mode", variable=self.dark_var, command=self._toggle_theme,
                        style="Switch.TCheckbutton" if sv_ttk else "TCheckbutton").pack(side="right", padx=6)

    # ---------------------------------------------------------------- theme
    def _reg(self, widget, **roles):
        """Record which palette colour each option of a widget uses (e.g. bg="surface", fg="muted"),
        so a theme change can repaint it. sv-ttk resets plain tk widgets to its own colours on every
        theme change, so colours are always re-applied from these roles afterwards."""
        self._themed.append((widget, roles))
        self._paint_one(widget, roles)
        return widget

    @staticmethod
    def _paint_one(widget, roles):
        pal = THEMES[THEME]
        try:
            widget.configure(**{opt: pal[key] for opt, key in roles.items()})
        except tk.TclError:  # widget destroyed
            pass

    def _paint(self):
        self.configure(bg=BG)
        for widget, roles in self._themed:
            self._paint_one(widget, roles)
        for text in self.texts:
            text.tag_configure("head", foreground=FG)
            text.tag_configure("warn", foreground=WARN)
        for tb in self.toolbars:
            self._style_toolbar(tb)

    def _toggle_theme(self):
        self._apply_theme("dark" if self.dark_var.get() else "light")

    def _apply_theme(self, name: str):
        _use_palette(name)
        self._style()
        self._paint()
        self.after_idle(self._paint)  # again once sv-ttk has finished repainting tk widgets
        for fig, _ in self.figs.values():
            fig.set_facecolor(SURFACE)
        if self.result:
            self._draw_all()
        else:
            for key in self.figs:
                self._placeholder(key)
        self._titlebar()
        _save_theme(name)

    def _set_status(self, text: str, role: str = "muted"):
        self.status.configure(text=text)
        for widget, roles in self._themed:
            if widget is self.status:
                roles["fg"] = role
                self._paint_one(widget, roles)

    def _style_toolbar(self, tb: NavigationToolbar2Tk):
        for w in (tb, *tb.winfo_children()):
            for opt, val in (("background", BG), ("foreground", FG), ("activebackground", BORDER),
                             ("selectcolor", BORDER)):
                try:
                    w.configure(**{opt: val})
                except tk.TclError:
                    pass
        try:  # re-render icons: matplotlib inverts them when the button background is dark
            for button in tb._buttons.values():
                tb._set_image_for_button(button)
        except (AttributeError, tk.TclError):
            pass

    def _titlebar(self):
        """Match the Windows title bar to the theme (Windows 10 20H1+ / 11)."""
        if os.name != "nt":
            return
        try:
            self.update_idletasks()
            hwnd = ctypes.windll.user32.GetParent(self.winfo_id())
            value = ctypes.c_int(1 if THEME == "dark" else 0)
            ctypes.windll.dwmapi.DwmSetWindowAttribute(hwnd, 20, ctypes.byref(value), ctypes.sizeof(value))
        except (AttributeError, OSError):
            pass

    # ---------------------------------------------------------------- input helpers
    def _section(self, parent, text):
        lf = ttk.LabelFrame(parent, text=text, padding=(12, 6, 12, 8))
        lf.pack(fill="x", padx=6, pady=6)
        lf.columnconfigure(0, weight=1)
        return lf

    def _field(self, parent, row, key, label, default, unit, tip):
        var = tk.StringVar(value=_fmt(default))
        lab = ttk.Label(parent, text=label)
        lab.grid(row=row, column=0, sticky="w", pady=2)
        ent = ttk.Entry(parent, textvariable=var, width=11, justify="right")
        ent.grid(row=row, column=1, sticky="e", padx=4)
        self._reg(ttk.Label(parent, text=unit, width=7), foreground="muted").grid(row=row, column=2, sticky="w")
        Tooltip(lab, tip)
        Tooltip(ent, tip)
        self._add_var(key, var, label, default)
        return ent

    def _combo(self, parent, row, key, label, values, default, tip, width=22):
        var = tk.StringVar(value=default)
        lab = ttk.Label(parent, text=label)
        lab.grid(row=row, column=0, sticky="w", pady=2)
        cb = ttk.Combobox(parent, textvariable=var, values=values, width=width, font=F)
        cb.grid(row=row, column=1, columnspan=2, sticky="e", padx=4)
        Tooltip(lab, tip)
        Tooltip(cb, tip)
        self._add_var(key, var, label, default)
        return cb

    def _check(self, parent, row, key, label, default, tip):
        var = tk.BooleanVar(value=default)
        cb = ttk.Checkbutton(parent, text=label, variable=var)
        cb.grid(row=row, column=0, columnspan=3, sticky="w", pady=2)
        Tooltip(cb, tip)
        self._add_var(key, var, label, default)
        return cb

    def _radios(self, parent, row, key, label, options, default, tip):
        var = tk.StringVar(value=default)
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=(2, 0))
        fr = ttk.Frame(parent)
        fr.grid(row=row + 1, column=0, columnspan=3, sticky="w", pady=(0, 2))
        for val, text in options:
            rb = ttk.Radiobutton(fr, text=text, value=val, variable=var)
            rb.pack(side="left", padx=(0, 10))
            Tooltip(rb, tip)
        self._add_var(key, var, label, default)

    def _note(self, parent, text, row=None):
        lbl = ttk.Label(parent, text=text, style="Muted.TLabel", wraplength=self.px(410), justify="left")
        if row is None:
            lbl.pack(anchor="w", padx=10, pady=(4, 8))
        else:
            lbl.grid(row=row, column=0, columnspan=3, sticky="w", pady=(2, 2))
        return lbl

    def _add_var(self, key, var, label, default):
        self.vars[key] = var
        self.labels[key] = label
        self.defaults[key] = _fmt(default) if not isinstance(default, bool) else default

    def _tab(self, nb, title):
        sf = ScrollFrame(nb)
        self._reg(sf.canvas, bg="bg")
        nb.add(sf, text=title)
        self.scrollers.append(sf)
        return sf.inner

    def _wheel(self, e):
        w = str(e.widget)
        if not w.startswith(str(self.inputs_nb)) or "combobox" in w.lower():
            return
        sf = self.scrollers[self.inputs_nb.index(self.inputs_nb.select())]
        sf.canvas.yview_scroll(-3 if e.delta > 0 else 3, "units")

    # ---------------------------------------------------------------- inputs
    def _build_inputs(self, parent):
        e, t, j = eng.EngineInputs(), rt.ThermalInputs(), inj.InjectorInputs()
        self.inputs_nb = nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)

        # ---- Engine
        tab = self._tab(nb, "Engine")
        s = self._section(tab, "Example cases")
        self._combo(s, 0, "preset", "Load", list(PRESETS), list(PRESETS)[0],
                    "Fill the inputs with an example case. 'Defaults' restores the dataclass defaults.",
                    width=28).bind("<<ComboboxSelected>>", self._apply_preset)

        s = self._section(tab, "Operating point")
        self._field(s, 0, "pc", "Chamber pressure", e.pc_bar, "bar",
                    "Stagnation pressure at the injector face.")
        self._field(s, 1, "thrust", "Thrust", e.thrust_n, "N",
                    "Thrust at the design ambient pressure.")
        self._field(s, 2, "pamb", "Ambient pressure", e.pamb_bar, "bar",
                    "Design ambient pressure. 1.01325 = sea level, 0 = vacuum (then give an expansion ratio).")
        qb = ttk.Frame(s)
        qb.grid(row=3, column=0, columnspan=3, sticky="w", pady=(0, 2))
        ttk.Button(qb, text="Sea level", width=10,
                   command=lambda: self.vars["pamb"].set("1.01325")).pack(side="left", padx=(0, 4))
        ttk.Button(qb, text="Vacuum", width=10, command=lambda: self.vars["pamb"].set("0")).pack(side="left")
        self._field(s, 4, "pe", "Exit pressure", e.pe_bar, "bar",
                    "Nozzle exit pressure. 'auto' = equal to ambient (ideal expansion).")
        self._field(s, 5, "eps", "Expansion ratio Ae/At", e.eps, "-",
                    "Exit area / throat area. 'auto' = solved from the exit pressure. A value here overrides it.")

        s = self._section(tab, "Propellants")
        self._combo(s, 0, "ox", "Oxidizer", OXIDIZERS, e.ox, "Oxidizer name as used by NASA CEA.")
        self._combo(s, 1, "fuel", "Fuel", FUELS, e.fuel, "Fuel name as used by NASA CEA. 'IPA' is liquid "
                    "isopropyl alcohol (this tool's card; RocketCEA's 'Isopropanol' card uses the gas-phase "
                    "heat of formation).")
        self._field(s, 2, "of", "Mixture ratio O/F (core)", e.of_ratio, "-",
                    "CORE oxidizer mass flow / CORE fuel mass flow. Film-cooling fuel is not included "
                    "(set the film on the Cooling tab, as a % of the core flow).")
        mode = "frozen_at_throat" if e.frozen_at_throat else ("frozen" if e.frozen else "equilibrium")
        self._radios(s, 3, "mode", "Nozzle chemistry", [("equilibrium", "Equilibrium"), ("frozen", "Frozen"),
                     ("frozen_at_throat", "Frozen at throat")], mode,
                     "Shifting equilibrium (upper bound), frozen from the chamber (lower bound), "
                     "or equilibrium to the throat and frozen after.")

        s = self._section(tab, "Mixture ratio sweep")
        sw = e.of_sweep or (max(e.of_ratio * 0.5, 0.5), e.of_ratio * 2.0, 0.1)
        self._check(s, 0, "sweep_on", "Run an O/F sweep", e.of_sweep is not None,
                    "Evaluate delivered performance over a range of O/F and plot it.")
        self._field(s, 1, "of_min", "O/F min", round(sw[0], 2), "-", "Start of the sweep.")
        self._field(s, 2, "of_max", "O/F max", round(sw[1], 2), "-", "End of the sweep.")
        self._field(s, 3, "of_step", "O/F step", sw[2], "-", "Sweep step size.")
        self._check(s, 4, "opt_of", "Design at the O/F with the highest delivered Isp", e.optimize_of,
                    "Overrides the O/F above with the best O/F from the sweep.")

        # ---- Chamber & nozzle
        tab = self._tab(nb, "Geometry")
        s = self._section(tab, "Chamber")
        self._field(s, 0, "lstar", "Characteristic length L*", e.lstar_m, "m",
                    "Chamber volume / throat area. Typically 0.8-1.5 m.")
        self._field(s, 1, "cr", "Contraction ratio Ac/At", e.contraction_ratio, "-",
                    "Chamber area / throat area. 'auto' = Humble's empirical correlation.")
        self._field(s, 2, "conv_angle", "Convergent half angle", e.conv_half_angle_deg, "deg",
                    "Half angle of the convergent cone. Typically 30-45 deg.")
        s = self._section(tab, "Nozzle (conical)")
        self._field(s, 0, "cone_angle", "Divergent half angle", e.cone_half_angle_deg, "deg",
                    "Half angle of the divergent cone. 15 deg is standard.")
        s = self._section(tab, "Efficiencies")
        self._field(s, 0, "eta_cstar", "c* efficiency", e.eta_cstar, "-",
                    "Combustion / injector efficiency applied to ideal c*. Typically 0.90-0.98.")
        self._field(s, 1, "eta_nozzle", "Nozzle efficiency", e.eta_nozzle, "-",
                    "Nozzle losses other than divergence (boundary layer, kinetics). Typically 0.97-0.99.")
        self._field(s, 2, "eta_div", "Divergence efficiency", e.eta_divergence, "-",
                    "'auto' = (1 + cos α) / 2 for the cone angle. A value here overrides it.")

        # ---- Cooling
        tab = self._tab(nb, "Cooling")
        s = self._section(tab, "Regenerative cooling")
        self._check(s, 0, "thermal_on", "Run the cooling analysis", True,
                    "Untick to size the engine and injector only.")
        s = self._section(tab, "Coolant (fuel)")
        cool_default = next((k for k, v in COOLANTS.items()
                             if v == rt.COOLANT_ALIASES.get(t.coolant, t.coolant)), t.coolant)
        self._combo(s, 0, "coolant", "Property model", list(COOLANTS), cool_default,
                    "The coolant is the fuel at the full fuel flow. 'Auto' selects the CoolProp fluid from the "
                    "fuel name (RP-1 uses the n-Dodecane surrogate; IPA uses the built-in isopropanol "
                    "correlations, as CoolProp has none). Any CoolProp fluid name can be typed.")
        self._field(s, 1, "t_in", "Inlet temperature", t.t_in_K, "K",
                    "Fuel temperature entering the cooling jacket.")
        self._check(s, 2, "counterflow", "Counterflow (enters at the nozzle exit)", t.counterflow,
                    "Ticked: fuel enters at the nozzle exit and flows toward the injector. "
                    "Unticked: enters at the injector end.")
        tf = ttk.Frame(s)
        tf.grid(row=3, column=0, columnspan=3, sticky="ew", pady=2)
        ttk.Label(tf, text="Property table").pack(side="left")
        self.vars["table"] = tk.StringVar(value=t.coolant_table or "")
        self.labels["table"] = "Property table"
        self.defaults["table"] = t.coolant_table or ""
        ent = ttk.Entry(tf, textvariable=self.vars["table"], width=20)
        ent.pack(side="left", padx=4)
        ttk.Button(tf, text="Browse…", command=self._browse_table).pack(side="left")
        ttk.Button(tf, text="Clear", command=lambda: self.vars["table"].set("")).pack(side="left", padx=(4, 0))
        Tooltip(ent, "Optional CSV with columns T_K, rho_kg_m3, mu_Pa_s, cp_J_kgK, k_W_mK. "
                     "Overrides the property model above.")
        self._note(s, "Inlet pressure is solved so the fuel reaches the injector at Pc × (1 + stiffness).", row=4)

        sd = ws.StressInputs()
        s = self._section(tab, "Liner wall")
        self._combo(s, 0, "material", "Material", list(MATERIALS), mat.get(sd.material).label,
                    "Liner material. Selecting one fills in its room-temperature conductivity below; the wall "
                    "stress check uses its temperature-dependent E, α, ν and strengths (see README for "
                    "sources).").bind("<<ComboboxSelected>>", self._apply_material)
        self._field(s, 1, "k_wall", "Thermal conductivity", t.k_wall, "W/m·K",
                    "Liner conductivity, constant along the wall (room-temperature value). "
                    "Copper alloys ~300-390, steels and nickel alloys ~11-25.")
        self._field(s, 2, "t_limit", "Hot-wall temperature limit", t.t_wall_limit_K, "K",
                    "Used for the margin check and warnings.")
        self._field(s, 3, "rough", "Channel roughness", t.roughness_um, "µm",
                    "Channel surface roughness for friction. Machined ~1-6 µm, additively made ~10-30 µm.")

        s = self._section(tab, "Wall stress")
        self._check(s, 0, "stress_on", "Check wall stresses and safety factors", True,
                    "Pressure bending between ribs, through-wall and longitudinal thermal stress (NASA SP-125), "
                    "checked with signs at the hot and cold face, rib root and mid-span.")
        self._field(s, 1, "sf_y", "Yield SF target (von Mises)", sd.sf_yield_target, "-",
                    "Minimum yield safety factor before the summary tile turns red.")
        self._field(s, 2, "t_jacket", "Jacket outer temperature", sd.t_jacket_K, "K",
                    "Temperature of the jacket (outer shell) that restrains the gas-side wall axially. "
                    "Longitudinal stress = E α (T_wg − T_jacket). Default: ambient.")
        self._field(s, 3, "life_factor", "Fatigue life factor", sd.life_factor, "-",
                    "Allowable firings = cycles to failure / this factor (scatter factor on mean-curve life).")
        self._note(s, "All structural results assume fully elastic behaviour, even past yield. Stresses above "
                      "yield are elastic-equivalent values; the buckling margin is optimistic where the wall "
                      "yields.", row=4)

        s = self._section(tab, "Model settings")
        self._field(s, 0, "curv", "Throat curvature R/Rt (Bartz)", t.throat_curv_ratio, "-",
                    "Throat radius of curvature used in Bartz's curvature term (the contour throat is sharp).")
        self._check(s, 1, "t0_eff", "Scale gas temperature by c* efficiency²", t.use_cstar_eff_T0,
                    "Ticked: T0 = Tc × η_c*². Unticked: ideal CEA chamber temperature.")
        self._radios(s, 2, "corr", "Coolant correlation", [("sieder-tate", "Sieder-Tate"),
                     ("dittus-boelter", "Dittus-Boelter")], t.correlation,
                     "Coolant-side Nusselt correlation. Sieder-Tate includes the wall viscosity correction.")
        self._field(s, 4, "n_stations", "Axial stations", t.n_stations, "-",
                    "Number of axial calculation stations.")

        s = self._section(tab, "Film cooling (fuel film from the injector face)")
        self._field(s, 0, "film_pct", "Film flow, % of CORE flow", e.film_pct_core, "%",
                    "Fuel film flow as a percentage of the CORE mass flow (core oxidizer + core fuel). "
                    "Example: a 1 kg/s core with 0.15 kg/s of film is 15 %. 0 = no film. The film is taken from "
                    "the fuel manifold after the regen jacket and adds mass flow but no thrust.")
        self._field(s, 1, "film_eff", "Film efficiency", t.film_efficiency, "-",
                    "Fraction of the film that stays on the wall (1 = all of it; no entrainment is modelled). "
                    "Lower it to anchor to test data.")
        self._field(s, 2, "film_h", "Gas-to-film h multiplier", t.film_h_scale, "× Bartz",
                    "Scales the gas-to-film heat-transfer coefficient. A recommended starting value is 1.8. Lucas & Golladay (NASA TN D-1988) needed a "
                    "somewhat higher value than 1 in the convergent section and throat.")
        self._note(s, "Liquid film (Grisson, non-entrained, non-reacting), then gaseous film (Hatch & Papell, "
                      "no mixing, no reaction). See README.", row=3)

        s = self._section(tab, "Anchoring multipliers")
        self._field(s, 0, "hg_scale", "Gas-side h_g multiplier", t.hg_scale, "× Bartz",
                    "Scales the Bartz gas-side coefficient everywhere, e.g. 0.7 to anchor to test data "
                    "or a CFD estimate. 1.0 = Bartz as published.  PDMS may apply a scalar of ~0.8.")
        self._field(s, 1, "hc_scale", "Coolant-side h_c multiplier", t.hc_scale, "×",
                    "Scales the coolant-side coefficient from the Nusselt correlation (applied before the "
                    "fin correction). 1.0 = correlation as published.")

        # ---- Channels
        tab = self._tab(nb, "Channels")
        tbl = sorted(t.channel_table)
        cyl, thr, ex = tbl[0][1:], min(tbl, key=lambda r_: abs(r_[0]))[1:], tbl[-1][1:]
        s = self._section(tab, "Channel count")
        self.n_ch_entry = self._field(s, 0, "n_ch", "Number of channels", t.n_channels, "-",
                                      "Channels around the circumference. Rib (land) width is what remains.")
        s = self._section(tab, "Auto-size widths")
        self.auto_w_cb = self._check(
            s, 0, "auto_w", "Size widths from the rib thicknesses", False,
            "Computes the channel width at the cylinder, throat and exit from the engine radii: "
            "w = 2π(r + wall t)/N − rib. Enter the rib at each section in the Rib column; height and "
            "wall thickness stay as entered.")
        self.auto_w_cb.configure(command=self._toggle_auto_width)
        self._note(s, "When ticked, the Rib column is editable and the widths are filled in after each run. "
                      "Ribs vary linearly between sections, like the other dimensions.", row=1)
        self.width_entries: list[ttk.Entry] = []
        self.height_entries: list[ttk.Entry] = []
        self.rib_entries: list[ttk.Entry] = []
        s = self._section(tab, "Channel dimensions (mm)")
        for c, h in enumerate(["", "Width", "Height", "Wall t", "Rib"]):
            ttk.Label(s, text=h, style="Head.TLabel").grid(row=0, column=c, padx=4, pady=(0, 2))
        for r_, (key, name, vals, tip) in enumerate([
                ("cyl", "Cylinder", (*cyl, 1.0), "Held constant along the chamber barrel."),
                ("thr", "Throat", (*thr, 1.0), "At the throat."),
                ("exit", "Nozzle exit", (*ex, 1.0), "At the nozzle exit.")], start=1):
            lab = ttk.Label(s, text=name)
            lab.grid(row=r_, column=0, sticky="w")
            Tooltip(lab, tip + " Interpolated linearly from the end of the cylinder to the throat, "
                               "then to the exit.")
            for c, (sub, v, what) in enumerate(zip(("w", "h", "t", "rib"), vals,
                                                   ("channel width", "channel height (rib height)",
                                                    "hot-wall thickness",
                                                    "rib (land) thickness, used when auto-sizing widths")),
                                               start=1):
                var = tk.StringVar(value=_fmt(v))
                en = ttk.Entry(s, textvariable=var, width=7, justify="right")
                en.grid(row=r_, column=c, padx=4, pady=2)
                Tooltip(en, f"{name}: {what} [mm]")
                self._add_var(f"{key}_{sub}", var, f"{name} {what}", v)
                if sub == "rib":
                    self.rib_entries.append(en)
                if sub == "w":
                    self.width_entries.append(en)
                if sub == "h":
                    self.height_entries.append(en)
        self._note(tab, "Cylinder dimensions are constant up to the start of convergence, then vary linearly "
                        "to the throat values and on to the exit values. Rib width = circumference / N − channel "
                        "width; results are listed under Summary. The run stops if the channels do not fit.")

        o = ro.OptimizerInputs()
        s = self._section(tab, "Channel optimizer")
        self._check(s, 0, "opt_on", "Optimize channels when running", False,
                    "Off: the channels above are used as entered. On: Run analysis first searches the number "
                    "of channels and the width and height at each section for the lowest peak hot-wall "
                    "temperature within the limits below, writes them into the table, then runs the analysis. "
                    "Wall thickness stays as entered.").configure(command=self._toggle_auto_width)
        self._note(s, "Takes several minutes. The current table is the starting point; progress is shown in "
                      "the status bar.", row=1)
        self.opt_entries = [
            self._field(s, 2, "opt_nmin", "Channels, minimum", o.n_min, "-", "Smallest channel count searched."),
            self._field(s, 3, "opt_nmax", "Channels, maximum", o.n_max, "-",
                        "Largest channel count searched. The search also stops where the minimum rib no longer fits."),
            self._field(s, 4, "opt_wmin", "Minimum width", o.w_min_mm, "mm", "Smallest channel width allowed."),
            self._field(s, 5, "opt_hmin", "Minimum height", o.h_min_mm, "mm", "Smallest channel height allowed."),
            self._field(s, 6, "opt_hmax", "Maximum height", o.h_max_mm, "mm", "Largest channel height allowed."),
            self._field(s, 7, "opt_rib", "Minimum rib", o.rib_min_mm, "mm",
                        "Smallest rib (land) width allowed anywhere along the jacket."),
            self._field(s, 8, "opt_dp", "Maximum regen ΔP", o.dp_max_bar, "bar",
                        "Largest allowed friction pressure drop through the cooling channels."),
        ]
        self._toggle_auto_width()

        # ---- Injector
        tab = self._tab(nb, "Injector")
        s = self._section(tab, "Injector stiffness (ΔP / Pc)")
        self._field(s, 0, "stiff_f", "Fuel", j.stiff_fuel, "-",
                    "Fuel injector pressure drop as a fraction of Pc. Typically 0.15-0.30.")
        self._field(s, 1, "stiff_o", "Oxidizer", j.stiff_ox, "-",
                    "Oxidizer injector pressure drop as a fraction of Pc. Typically 0.15-0.30.")
        s = self._section(tab, "Orifices (plain holes)")
        for c, h in enumerate(["", "Fuel", "Oxidizer", "Film"]):
            ttk.Label(s, text=h, style="Head.TLabel").grid(row=0, column=c, padx=4, pady=(0, 2))
        side_names = {"fuel": "Fuel", "ox": "Oxidizer", "film": "Film-cooling holes (fuel manifold; used when film > 0)"}
        for r_, (key, label, vals, tip) in enumerate([
                ("n", "Number of holes", (j.n_fuel, j.n_ox, j.n_film), "Number of identical orifices on each side."),
                ("cd", "Discharge coefficient", (j.cd_fuel, j.cd_ox, j.cd_film),
                 "Sharp-edged drilled hole ~0.60-0.65, chamfered or rounded inlet ~0.8.")], start=1):
            lab = ttk.Label(s, text=label)
            lab.grid(row=r_, column=0, sticky="w")
            Tooltip(lab, tip)
            for c, (side, v) in enumerate(zip(("fuel", "ox", "film"), vals), start=1):
                var = tk.StringVar(value=_fmt(v))
                en = ttk.Entry(s, textvariable=var, width=9, justify="right")
                en.grid(row=r_, column=c, padx=4, pady=2)
                Tooltip(en, f"{side_names[side]}: {tip}")
                self._add_var(f"{key}_{side}", var, f"{side} {label.lower()}", v)
        s = self._section(tab, "Oxidizer at the injector")
        self._field(s, 0, "t_ox", "Temperature", j.t_ox_K, "K",
                    "Oxidizer temperature in the injector manifold. LOX ~90 K, N2O ~280 K.")
        self._field(s, 1, "rho_ox", "Density", j.ox_density, "kg/m³",
                    "'auto' = CoolProp from the oxidizer name (LOX, GOX, N2O, F2). Required for other oxidizers.")
        self._note(tab, "Fuel path: regen channels → injector manifold → orifices → chamber. "
                        "Oxidizer path: manifold → orifices → chamber. Manifold pressure = Pc × (1 + stiffness). "
                        "Fuel density at the orifices is taken at the regen outlet. "
                        "Feed lines and pumps upstream are not included.")

    def _build_buttons(self, parent):
        fr = ttk.Frame(parent)
        fr.pack(fill="x", pady=(8, 0))
        self.go_btn = ttk.Button(fr, text="Run analysis", style="Accent.TButton" if sv_ttk else "Run.TButton",
                                 cursor="hand2", command=self.run)
        self.go_btn.pack(fill="x", padx=6, ipady=4)
        row = ttk.Frame(fr)
        row.pack(fill="x", pady=6, padx=6)
        ttk.Button(row, text="Export results…", command=self._export).pack(side="left")
        ttk.Button(row, text="Reset to defaults", command=self._reset).pack(side="left", padx=6)
        row = ttk.Frame(fr)
        row.pack(fill="x", pady=(0, 6), padx=6)
        ttk.Button(row, text="Save settings…", command=self._save_settings).pack(side="left")
        ttk.Button(row, text="Load settings…", command=self._load_settings).pack(side="left", padx=6)

    # ---------------------------------------------------------------- outputs
    def _build_outputs(self, parent):
        self.out_nb = nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)

        stats = ttk.Frame(nb, padding=(6, 6, 6, 0))
        nb.add(stats, text="Summary")
        self.cards = {}
        self._cards(stats, CARD_NAMES)
        self.text = self._textbox(stats)
        self._say("Set the inputs on the left, or load an example case, and click Run analysis.", tag="head")

        # Structural results: their own section (tiles + report), and their own section of plots
        structs = ttk.Frame(nb, padding=(6, 6, 6, 0))
        nb.add(structs, text="Structures")
        self._reg(ttk.Label(structs, text=ELASTIC_NOTE, wraplength=self.px(900), justify="left"),
                  foreground="warn").pack(anchor="w", padx=4, pady=(0, 6))
        self._cards(structs, STRUCT_CARD_NAMES)
        self.stext = self._textbox(structs)
        self._say("Wall stresses, safety factors, fatigue life and buckling appear here after a run "
                  "(Cooling tab → Wall stress).", tag="head", widget=self.stext)
        splots = ttk.Frame(nb, padding=(0, 4, 0, 0))
        nb.add(splots, text="Structural plots")
        self._reg(ttk.Label(splots, text=ELASTIC_NOTE, wraplength=self.px(900), justify="left"),
                  foreground="warn").pack(anchor="w", padx=10, pady=(2, 4))
        self.struct_nb = ttk.Notebook(splots)
        self.struct_nb.pack(fill="both", expand=True)

        for key, title, book in [("engine", "Engine", nb), ("temps", "Temperatures", nb), ("flux", "Heat flux", nb),
                                 ("coolant", "Coolant", nb),
                                 ("stress", "Envelope & yield SF", self.struct_nb),
                                 ("components", "Stress components", self.struct_nb),
                                 ("life", "Fatigue & buckling", self.struct_nb),
                                 ("budget", "Pressure budget", nb), ("sweep", "O/F sweep", nb)]:
            fr = ttk.Frame(book)
            book.add(fr, text=title)
            fig = Figure(figsize=(9, 6), dpi=96 * self.ui_scale, layout="tight")
            canvas = FigureCanvasTkAgg(fig, master=fr)
            tb = NavigationToolbar2Tk(canvas, fr, pack_toolbar=False)
            tb.update()
            self._style_toolbar(tb)
            self.toolbars.append(tb)
            tb.pack(side="bottom", fill="x")
            canvas.get_tk_widget().configure(highlightthickness=0)
            self._reg(canvas.get_tk_widget(), bg="surface")
            canvas.get_tk_widget().pack(fill="both", expand=True)
            self.figs[key] = (fig, canvas)
            # A plot is visible when its own tab is selected and, if nested, the outer tab holding it too
            self.fig_tabs[key] = (book, str(fr), None if book is nb else str(splots))
            self._placeholder(key)
        for book in (nb, self.struct_nb):
            book.bind("<<NotebookTabChanged>>", lambda e: self._draw_stale())

    def _cards(self, parent, names):
        """A grid of summary tiles (6 per row) registered in self.cards."""
        cards = ttk.Frame(parent)
        cards.pack(fill="x", pady=(0, 6))
        for i, name in enumerate(names):
            c = self._reg(tk.Frame(cards, highlightthickness=1, padx=14, pady=10),
                          bg="surface", highlightbackground="border")
            c.grid(row=i // 6, column=i % 6, sticky="nsew", padx=3, pady=3)
            cards.columnconfigure(i % 6, weight=1, uniform="cards")
            self._reg(tk.Label(c, text=name, font=(UI, 9), anchor="w"), bg="surface", fg="muted").pack(fill="x")
            val_roles = {"bg": "surface", "fg": "fg"}
            val = self._reg(tk.Label(c, text="—", font=(UI_DISPLAY, 15, "bold"), anchor="w"), **val_roles)
            val_roles = self._themed[-1][1]  # keep the stored dict so _card can change the colour role
            val.pack(fill="x")
            sub = self._reg(tk.Label(c, text="", font=(UI, 9), anchor="w"), bg="surface", fg="muted")
            sub.pack(fill="x")
            self.cards[name] = (val, sub, val_roles)

    def _textbox(self, parent) -> tk.Text:
        tf = ttk.Frame(parent)
        tf.pack(fill="both", expand=True)
        text = self._reg(tk.Text(tf, font=MONO, wrap="none", bd=0, highlightthickness=1, padx=8, pady=6),
                         bg="surface", fg="fg", highlightbackground="border", insertbackground="fg")
        ys = ttk.Scrollbar(tf, command=text.yview)
        xs = ttk.Scrollbar(parent, orient="horizontal", command=text.xview)
        text.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)
        ys.pack(side="right", fill="y")
        text.pack(side="left", fill="both", expand=True)
        xs.pack(fill="x")
        text.tag_configure("head", font=(UI_DISPLAY, 12, "bold"), foreground=FG, spacing1=8, spacing3=4)
        text.tag_configure("warn", foreground=WARN)
        self.texts.append(text)
        return text

    def _placeholder(self, key, msg="Run the analysis to see results."):
        fig, canvas = self.figs[key]
        fig.clear()
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11, color=MUTED, transform=ax.transAxes)
        ax.axis("off")
        self._render(key)

    def _render(self, key):
        """Draw a plot now if its tab is showing, otherwise when the tab is next selected. A hidden
        canvas has no real size yet, so drawing it early is wasted and is redone at the right size."""
        self.stale.add(key)
        self._draw_stale()

    def _draw_stale(self):
        outer = self.out_nb.select()
        for key in list(self.stale):
            book, frame, holder = self.fig_tabs[key]
            if book.select() == frame and (holder is None or outer == holder):
                self.stale.discard(key)
                self.figs[key][1].draw_idle()

    def _say(self, msg, tag=None, clear=True, widget=None):
        text = widget or self.text
        text.configure(state="normal")
        if clear:
            text.delete("1.0", "end")
        text.insert("end", msg, tag)
        text.configure(state="disabled")

    # ---------------------------------------------------------------- input parsing
    def _num(self, key, optional=False, positive=True) -> float | None:
        s = str(self.vars[key].get()).strip()
        if s == "" or s.lower() in ("auto", "none"):
            if optional:
                return None
            raise ValueError(f"'{self.labels[key]}' is required.")
        try:
            v = float(s)
        except ValueError:
            raise ValueError(f"'{self.labels[key]}' must be a number (got '{s}').") from None
        if positive and v <= 0:
            raise ValueError(f"'{self.labels[key]}' must be greater than zero (got {v:g}).")
        return v

    def _collect(self):
        g, v = self._num, self.vars
        sweep = None
        if v["sweep_on"].get():
            sweep = (g("of_min"), g("of_max"), g("of_step"))
            if sweep[1] <= sweep[0]:
                raise ValueError("O/F sweep: the maximum must be greater than the minimum.")
        for k in ("eta_cstar", "eta_nozzle"):
            if g(k) > 1:
                raise ValueError(f"'{self.labels[k]}' cannot exceed 1.")
        mode = v["mode"].get()
        ei = eng.EngineInputs(
            pc_bar=g("pc"), thrust_n=g("thrust"), pamb_bar=g("pamb", positive=False),
            pe_bar=g("pe", True), eps=g("eps", True),
            ox=v["ox"].get().strip(), fuel=v["fuel"].get().strip(), of_ratio=g("of"),
            film_pct_core=g("film_pct", positive=False), of_sweep=sweep, optimize_of=bool(sweep) and v["opt_of"].get(),
            frozen=mode == "frozen", frozen_at_throat=mode == "frozen_at_throat",
            lstar_m=g("lstar"), contraction_ratio=g("cr", True), conv_half_angle_deg=g("conv_angle"),
            cone_half_angle_deg=g("cone_angle"), eta_cstar=g("eta_cstar"), eta_nozzle=g("eta_nozzle"),
            eta_divergence=g("eta_div", True), csv_path=None, plot=False, plot_path=None)
        if ei.pamb_bar < 0:
            raise ValueError("Ambient pressure cannot be negative.")
        if ei.film_pct_core < 0:
            raise ValueError("Film flow (% of core flow) cannot be negative.")

        cooling = v["thermal_on"].get()
        cool = v["coolant"].get().strip()
        th = rt.ThermalInputs(
            coolant=COOLANTS.get(cool, cool), t_in_K=g("t_in"), counterflow=v["counterflow"].get(),
            coolant_table=v["table"].get().strip() or None,
            n_channels=int(g("n_ch")), k_wall=g("k_wall"), roughness_um=g("rough", positive=False),
            t_wall_limit_K=g("t_limit"), throat_curv_ratio=g("curv"), use_cstar_eff_T0=v["t0_eff"].get(),
            hg_scale=g("hg_scale"), hc_scale=g("hc_scale"),
            film_efficiency=g("film_eff"), film_h_scale=g("film_h"),
            correlation=v["corr"].get(), n_stations=int(g("n_stations")),
            csv_path=None, plot=False, plot_path=None)
        opt = None
        if v["opt_on"].get():
            if not cooling:
                raise ValueError("The channel optimizer needs the cooling analysis switched on.")
            opt = ro.OptimizerInputs(
                n_min=int(g("opt_nmin")), n_max=int(g("opt_nmax")), w_min_mm=g("opt_wmin"),
                h_min_mm=g("opt_hmin"), h_max_mm=g("opt_hmax"), rib_min_mm=g("opt_rib"), dp_max_bar=g("opt_dp"))
        rib = tuple(g(f"{s}_rib") for s in ("cyl", "thr", "exit")) if v["auto_w"].get() and not opt else None
        sections = tuple((None if rib else g(f"{s}_w"), g(f"{s}_h"), g(f"{s}_t")) for s in ("cyl", "thr", "exit"))
        ji = inj.InjectorInputs(
            stiff_fuel=g("stiff_f"), stiff_ox=g("stiff_o"), n_fuel=int(g("n_fuel")), n_ox=int(g("n_ox")),
            cd_fuel=g("cd_fuel"), cd_ox=g("cd_ox"), t_ox_K=g("t_ox"), ox_density=g("rho_ox", True),
            n_film=int(g("n_film")), cd_film=g("cd_film"))
        for k in ("cd_fuel", "cd_ox", "cd_film", "film_eff"):
            if g(k) > 1:
                raise ValueError(f"'{self.labels[k]}' cannot exceed 1.")
        si = None
        if cooling and v["stress_on"].get():
            label = v["material"].get().strip()
            if label not in mat.BY_LABEL:
                raise ValueError("Select a liner material on the Cooling tab for the wall stress check "
                                 "(or untick 'Check wall stresses').")
            si = ws.StressInputs(material=mat.BY_LABEL[label].name, sf_yield_target=g("sf_y"),
                                 t_jacket_K=g("t_jacket"),
                                 life_factor=g("life_factor"))
        return ei, th, sections, ji, cooling, rib, opt, si

    # ---------------------------------------------------------------- run
    def run(self):
        if self.busy:
            return
        try:
            args = self._collect()
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return
        self.busy = True
        self.go_btn.configure(state="disabled", text="Running…")
        self._set_status("Optimizing channels…" if args[6] else "Running analysis…")
        threading.Thread(target=self._work, args=args, daemon=True).start()
        self.after(100, self._poll)

    def _work(self, ei, th, sections, ji, cooling, rib, opt, si):
        def capture(fn, *a):
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                fn(*a)
            return buf.getvalue()

        try:
            design = eng.design_engine(ei)
            widths, optimized = None, None
            if opt is not None:
                start = tuple(s[0] for s in sections) + tuple(s[1] for s in sections)
                optimized = ro.optimize_channels(design, th, ji, tuple(s[2] for s in sections), opt,
                                                 start_mm=start, progress=lambda m: self.jobs.put(("progress", m)))
                th, feed = optimized["th"], optimized["feed"]
            else:
                if rib is not None:
                    widths = rt.widths_for_ribs(design, th.n_channels, rib, tuple(s[2] for s in sections))
                    sections = tuple((w, h, t) for w, (_, h, t) in zip(widths, sections))
                th.channel_table = rt.channel_table_from_sections(design, *sections)
                feed = inj.solve_feed_system(design, th, ji, cooling=cooling)
            res = feed["thermal"]
            if res and si is not None:
                feed["stress"] = ws.wall_stresses(design, th, res, si)
            texts = dict(engine=capture(eng.print_report, design),
                         thermal=capture(rt.print_thermal_report, design, th, res) if res else "",
                         feed=capture(inj.print_feed_report, design, feed))
            if feed.get("stress"):
                texts["stress"] = capture(ws.print_stress_report, design, res, feed["stress"])
            if optimized:
                texts["optimizer"] = capture(ro.print_optimize_report, optimized)
            self.jobs.put(("ok", dict(design=design, th=th, res=res, feed=feed, texts=texts, widths=widths,
                                      optimized=optimized)))
        except Exception as exc:  # anything from CEA / CoolProp / geometry
            self.jobs.put(("err", exc))

    def _poll(self):
        try:
            kind, payload = self.jobs.get_nowait()
        except queue.Empty:
            self.after(100, self._poll)
            return
        if kind == "progress":
            self._set_status(payload)
            self.after(100, self._poll)
            return
        self.busy = False
        self.go_btn.configure(state="normal", text="Run analysis")
        if kind == "err":
            self._set_status("Analysis failed. Check the inputs.", "bad")
            messagebox.showerror("Analysis failed", f"{type(payload).__name__}: {payload}")
            return
        self.result = payload
        design, th, res = payload["design"], payload["th"], payload["res"]
        if payload["widths"]:
            for key, w in zip(("cyl", "thr", "exit"), payload["widths"]):
                self.vars[f"{key}_w"].set(f"{w:.3f}")
        if payload["optimized"]:
            self.vars["n_ch"].set(str(payload["optimized"]["n"]))
            for key, (w, h, _) in zip(("cyl", "thr", "exit"), payload["optimized"]["sections"]):
                self.vars[f"{key}_w"].set(f"{w:.3f}")
                self.vars[f"{key}_h"].set(f"{h:.3f}")
        self._fill_stats(design, th, res, payload["feed"], payload["texts"])
        self._draw_all()
        stress = payload["feed"].get("stress")
        n_warn = (len(payload["feed"]["warnings"]) + (len(res["warnings"]) if res else 0)
                  + (len(stress["warnings"]) if stress else 0))
        film_on = design.sizing["mdot_film"] > 0
        msg = f"Analysis complete. Delivered Isp {design.perf['isp_film' if film_on else 'isp']:.1f} s"
        if film_on:
            msg += " (with film)"
        if res:
            tmax = self._tmax(res)
            margin = th.t_wall_limit_K - tmax
            msg += f", max hot-wall T {tmax:.0f} K (" + (f"{margin:.0f} K below limit)" if margin >= 0
                                                          else f"{-margin:.0f} K over limit)")
        where = "Summary and Structures tabs" if stress else "summary"
        msg += f". {n_warn} warning{'s' if n_warn != 1 else ''} listed in the {where}." if n_warn else "."
        self._set_status(msg, "fg")

    @staticmethod
    def _tmax(res) -> float:
        return float(np.nanmax(np.where(res["cooled"], res["out"]["t_wg"], np.nan)))

    # ---------------------------------------------------------------- summary
    def _fill_stats(self, design, th, res, feed, texts):
        p, s, i = design.perf, design.sizing, design.inputs
        if s["mdot_film"] > 0:   # film adds mass, no thrust: delivered Isp = F / (core + film)
            self._card("Isp (design)", f"{p['isp_film']:.1f} s", f"with film; core {p['isp']:.1f} s")
            self._card("Isp (vacuum)", f"{p['isp_vac_film']:.1f} s", f"with film; core {p['isp_vac']:.1f} s")
        else:
            self._card("Isp (design)", f"{p['isp']:.1f} s", f"ideal {p['isp_ideal']:.1f} s")
            self._card("Isp (vacuum)", f"{p['isp_vac']:.1f} s", f"sea level {p['isp_sl']:.1f} s")
        self._card("Thrust", f"{i.thrust_n:.0f} N", f"vacuum {s['thrust_vac']:.0f} N")
        self._card("Mass flow", f"{s['mdot']:.3f} kg/s",
                   f"core O/F {p['of']:.2f}" + (f"; + film {s['mdot_film']:.3f} kg/s" if s["mdot_film"] > 0 else ""))
        self._card("c*", f"{p['cstar']:.0f} m/s", f"Tc {p['tc']:.0f} K")
        self._card("Throat diameter", f"{s['dt'] * 1e3:.2f} mm", f"chamber {s['dc'] * 1e3:.1f} mm")
        self._card("Exit diameter", f"{s['de'] * 1e3:.2f} mm", f"ε = {p['eps']:.2f}")
        if res:
            tmax = self._tmax(res)
            margin = th.t_wall_limit_K - tmax
            self._card("Max hot-wall T", f"{tmax:.0f} K",
                       f"{margin:.0f} K below limit" if margin >= 0 else f"{-margin:.0f} K over limit",
                       role="ok" if margin >= 0.1 * th.t_wall_limit_K else ("warn" if margin >= 0 else "bad"))
            self._card("Regen ΔP", f"{res['dp_regen'] / eng.BAR:.2f} bar",
                       f"fuel ΔT {res['t_out'] - th.t_in_K:.0f} K, Q {res['q_total'] / 1e3:.0f} kW")
        else:
            self._card("Max hot-wall T", "—", "cooling analysis off")
            self._card("Regen ΔP", "—", "cooling analysis off")
        b, fu, ox = feed["budget"], feed["fuel"], feed["ox"]
        self._card("Feed pressure", f"{b['p_feed_f'] / eng.BAR:.1f} / {b['p_feed_o'] / eng.BAR:.1f} bar",
                   "fuel (regen inlet) / oxidizer")
        self._card("Orifice diameter", f"{fu['d'] * 1e3:.2f} / {ox['d'] * 1e3:.2f} mm",
                   f"{fu['n']} fuel / {ox['n']} oxidizer holes")
        film = res.get("film") if res else None
        if s["mdot_film"] <= 0:
            self._card("Liquid film length", "—", "no film cooling")
        elif not film:
            self._card("Liquid film length", "—", f"film {i.film_pct_core:g} % of core; cooling analysis off")
        else:
            value = "whole engine" if film["liquid_to_exit"] else f"{film['l_c'] * 1e3:.1f} mm"
            fh = feed.get("film")
            holes = f"; {fh['n']} × {fh['d'] * 1e3:.2f} mm holes" if fh else ""
            self._card("Liquid film length", value, f"film {i.film_pct_core:g} % of core{holes}")
        self._fill_structures(feed.get("stress"), texts)

        self._say("Engine and performance\n", "head")
        self._say(texts["engine"], clear=False)
        if res:
            self._say("\nChannels and ribs\n", "head", clear=False)
            self._say(self._ribs_text(design, th, res), clear=False)
            self._say("\nThermal\n", "head", clear=False)
            self._say_with_warnings(texts["thermal"])
        if "stress" in texts:
            self._say("\nStructures\n", "head", clear=False)
            self._say("  Wall stresses, safety factors, fatigue life and buckling: see the Structures tab.\n",
                      clear=False)
        self._say("\nFeed system and injector\n", "head", clear=False)
        self._say_with_warnings(texts["feed"])
        if "optimizer" in texts:
            self._say("\nChannel optimizer\n", "head", clear=False)
            self._say_with_warnings(texts["optimizer"])

    def _fill_structures(self, sres, texts):
        """Structures tab: governing-value tiles and the full structural report."""
        if not sres:
            for name in STRUCT_CARD_NAMES:
                self._card(name, "—", "stress check off")
            self._say("Wall stress check is off (Cooling tab → Wall stress), or cooling analysis is off.",
                      tag="head", widget=self.stext)
            return
        ex, si, mpa = sres["extremes"], sres["inputs"], 1e-6

        def where(e):
            return f"{e['face']}, {e['position']}, {e['region']}, x {e['x'] * 1e3:.0f} mm"

        e = ex.get("vm")
        self._card("Max von Mises stress", f"{e['value'] * mpa:.0f} MPa" if e else "—", where(e) if e else "")
        e = ex.get("sp")
        self._card("Max principal stress", f"{e['value'] * mpa:.0f} MPa" if e else "—",
                   f"{e['component']}, {e['sense']}, {e['face']}, x {e['x'] * 1e3:.0f} mm" if e else "")
        e = ex.get("sfy")
        if e is None:
            self._card("Min yield SF", "—", "no yield-strength data")
        else:
            goal = si.sf_yield_target
            self._card("Min yield SF", f"{e['value']:.2f}", f"vs von Mises; {where(e)}",
                       role="ok" if e["value"] >= goal else ("warn" if e["value"] >= 1 else "bad"))
        e = ex.get("nf")
        if e:
            allow = e["value"] / si.life_factor
            self._card("Fatigue life", f"{e['value']:.3g} cycles",
                       f"{allow:.3g} firings allowed (÷{si.life_factor:g}); strain {e['eps'] * 100:.2f} %",
                       role="ok" if allow >= 100 else ("warn" if allow >= 10 else "bad"))
        else:
            self._card("Fatigue life", "—", "no fatigue data for this material")
        e = ex.get("bm")
        if e:
            self._card("Buckling margin", f"{e['value']:.2f}",
                       f"elastic (upper bound), {e['region']}, x {e['x'] * 1e3:.0f} mm",
                       role="ok" if e["value"] >= 0 else "bad")
        else:
            self._card("Buckling margin", "—", "")
        self._say("Structural analysis of the hot wall (fully elastic)\n", "head", widget=self.stext)
        self._say_with_warnings(texts.get("stress", ""), widget=self.stext)

    def _say_with_warnings(self, text, widget=None):
        body, _, warns = text.partition(" WARNINGS")
        self._say(body, clear=False, widget=widget)
        if warns:
            self._say(" WARNINGS" + warns, "warn", clear=False, widget=widget)

    def _card(self, name, value, sub, role="fg"):
        """Set a summary tile. role picks the value colour from the palette: fg, ok, warn or bad."""
        val, lab, roles = self.cards[name]
        roles["fg"] = role
        val.configure(text=value)
        self._paint_one(val, roles)
        lab.configure(text=sub)

    def _ribs_text(self, design, th, res) -> str:
        st, o = res["st"], res["out"]
        x = st["x"]
        idx = [("Cylinder", int(np.argmin(abs(x + design.chamber["l_conv"])))),
               ("Throat", int(np.argmin(abs(x)))), ("Exit", len(x) - 1)]
        n = th.n_channels
        lines = [f"  Channels: {n}   flow area / channel and total at each section\n",
                 f"  {'section':<10}{'width':>7}{'height':>8}{'wall t':>8}{'rib':>7}{'pitch':>7}{'h/w':>6}"
                 f"{'Dh':>7}{'A_tot':>8}{'V':>7}{'Re':>8}{'eta_fin':>8}{'h_c':>8}",
                 f"  {'':<10}{'mm':>7}{'mm':>8}{'mm':>8}{'mm':>7}{'mm':>7}{'':>6}{'mm':>7}{'mm2':>8}"
                 f"{'m/s':>7}{'':>8}{'':>8}{'kW/m2K':>8}"]
        for name, j in idx:
            w, h, t = res["w"][j] * 1e3, res["h"][j] * 1e3, res["t"][j] * 1e3
            land = o["land"][j] * 1e3
            lines.append(f"  {name:<10}{w:7.2f}{h:8.2f}{t:8.2f}{land:7.2f}{w + land:7.2f}{h / w:6.2f}"
                         f"{o['dh'][j] * 1e3:7.2f}{n * w * h:8.1f}{o['vel'][j]:7.2f}{o['re'][j]:8.0f}"
                         f"{o['eta_f'][j]:8.3f}{o['hc'][j] / 1e3:8.2f}")
        lines.append(f"\n  Thinnest rib: {np.nanmin(o['land']) * 1e3:.2f} mm   "
                     f"(rib height = channel height; ribs are treated as fins)\n")
        return "\n".join(lines) + "\n"

    # ---------------------------------------------------------------- plots
    def _draw_all(self):
        design, th, res, feed = (self.result[k] for k in ("design", "th", "res", "feed"))
        self._draw_engine(design, res)
        if res:
            self._draw_temps(th, res)
            self._draw_flux(res)
            self._draw_coolant(res)
        else:
            for k in ("temps", "flux", "coolant"):
                self._placeholder(k, "Cooling analysis is turned off.")
        if feed.get("stress"):
            self._draw_stress(res, feed["stress"])
            self._draw_components(res, feed["stress"])
            self._draw_life(res, feed["stress"])
        else:
            for k in ("stress", "components", "life"):
                self._placeholder(k, "Wall stress check is off (Cooling tab).")
        self._draw_budget(feed)
        self._draw_sweep(design)

    def _ax(self, key, title, nrows=1):
        fig, _ = self.figs[key]
        fig.clear()
        axes = fig.subplots(nrows, 1, sharex=True) if nrows > 1 else fig.add_subplot(111)
        (axes[0] if nrows > 1 else axes).set_title(title, loc="left")
        return fig, axes

    def _finish(self, key):
        self._render(key)  # the figure's tight layout engine sets the margins at draw time

    def _draw_engine(self, design, res):
        i, p, s = design.inputs, design.perf, design.sizing
        fig, ax = self._ax("engine", f"{i.ox}/{i.fuel}   Pc {i.pc_bar:g} bar   F {i.thrust_n:g} N   "
                                     f"Isp {p['isp']:.0f} s   Dt {s['dt'] * 1e3:.1f} mm")
        if res:
            x, r = res["st"]["x"] * 1e3, res["st"]["r"] * 1e3
            t, h, c = res["t"] * 1e3, res["h"] * 1e3, res["cooled"]
            ax.fill_between(x, -r, r, color=PC["gas"], label="hot gas")
            for sg in (1, -1):
                ax.fill_between(x, sg * r, sg * (r + t), color=PC["liner"], label="liner" if sg > 0 else None)
                ax.fill_between(x, sg * (r + t), sg * (r + t + h), where=c, color=PC["channel"],
                                label="coolant channel" if sg > 0 else None)
            ax.plot(x, r, color=FG, lw=1.2)
            ax.plot(x, -r, color=FG, lw=1.2)
        else:
            cx, cr = design.contour
            ax.fill_between(cx * 1e3, -cr * 1e3, cr * 1e3, color=PC["gas"], label="hot gas")
            ax.plot(cx * 1e3, cr * 1e3, cx * 1e3, -cr * 1e3, color=FG, lw=1.2)
        ax.axvline(0, color=BAD, ls="--", lw=0.9, label="throat")
        ax.axvline(design.chamber["x_inj"] * 1e3, color=ACCENT, ls="--", lw=0.9, label="injector face")
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]  (0 = throat)")
        ax.set_ylabel("r [mm]")
        ax.legend(loc="lower left")
        self._finish("engine")

    def _draw_temps(self, th, res):
        x, o = res["st"]["x"] * 1e3, res["out"]
        fig, ax = self._ax("temps", "Wall and coolant temperatures")
        ax.plot(x, o["t_wg"], color=PC["hot"], lw=1.8, label="hot wall T_wg")
        ax.plot(x, o["t_wc"], color=PC["warm"], lw=1.8, label="cold wall T_wc")
        ax.plot(x, o["t_b"], color=PC["cool"], lw=1.8, label="coolant T_b")
        ax.axhline(th.t_wall_limit_K, color=MUTED, ls=":", lw=1.2, label="hot-wall limit")
        ax2 = ax.twinx()
        ax2.plot(x, o["t_aw"], color=PC["aw"], ls="--", lw=1.1, label="gas recovery T_aw")
        ax2.set_ylabel("T_aw [K]", color=MUTED)
        ax2.grid(False)
        ax2.spines["right"].set_visible(True)
        j = int(np.nanargmax(np.where(res["cooled"], o["t_wg"], np.nan)))
        ax.plot([x[j]], [o["t_wg"][j]], "o", color=PC["hot"], ms=5)
        ax.annotate(f"max {o['t_wg'][j]:.0f} K", xy=(x[j], o["t_wg"][j]), xytext=(8, 6),
                    textcoords="offset points", fontsize=8, color=PC["hot"])
        ax.set_xlabel("x [mm]  (0 = throat)")
        ax.set_ylabel("T [K]")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="center left")
        self._finish("temps")

    def _draw_flux(self, res):
        x, o = res["st"]["x"] * 1e3, res["out"]
        fig, ax = self._ax("flux", "Heat flux and heat-transfer coefficients")
        ax.plot(x, o["q"] / 1e6, color=PC["violet"], lw=1.8, label="heat flux q''")
        ax.set_ylabel("q'' [MW/m²]")
        ax2 = ax.twinx()
        ax2.plot(x, o["hg"] / 1e3, color=PC["green"], ls="--", lw=1.3, label="gas side h_g (Bartz)")
        ax2.plot(x, o["hc"] / 1e3, color=PC["sky"], ls="-.", lw=1.3, label="coolant side h_c")
        ax2.set_ylabel("h [kW/m²K]")
        ax2.grid(False)
        ax2.spines["right"].set_visible(True)
        ax.set_xlabel("x [mm]  (0 = throat)")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="upper left")
        self._finish("flux")

    def _draw_coolant(self, res):
        x, o = res["st"]["x"] * 1e3, res["out"]
        fig, (a1, a2, a3) = self._ax("coolant", "Coolant side", nrows=3)
        a1.plot(x, o["vel"], color=PC["cool"], lw=1.6)
        a1.set_ylabel("velocity [m/s]")
        a2.plot(x, o["re"], color=PC["violet"], lw=1.6)
        a2.axhline(1e4, color=MUTED, ls=":", lw=1, label="Re = 10⁴ (turbulent correlation lower limit)")
        a2.set_ylabel("Reynolds number")
        a2.legend(loc="upper left")
        a3.plot(x, o["p_b"] / eng.BAR, color=PC["green"], lw=1.6)
        a3.set_ylabel("pressure [bar]")
        a3.set_xlabel("x [mm]  (0 = throat)")
        self._finish("coolant")

    def _draw_stress(self, res, sres):
        x, ex, mpa = res["st"]["x"] * 1e3, sres["extremes"], 1e-6
        fig, (a1, a2, a3) = self._ax("stress", f"Stress envelope and yield safety factor ({sres['material'].name}; "
                                               "worst of the 4 check points; fully elastic)", nrows=3)
        si = sres["inputs"]
        panels = [(a1, sres["vm"] * mpa, "von Mises σ_vm\n[MPa]", PC["hot"], "vm", None),
                  (a2, sres["sp"] * mpa, "max principal σ_p\n[MPa]", PC["violet"], "sp", None),
                  (a3, sres["sfy"], "yield SF\n(σ_y / σ_vm)", PC["green"], "sfy", si.sf_yield_target)]
        for ax, y, label, color, key, target in panels:
            ax.plot(x, y, color=color, lw=1.6)
            ax.set_ylabel(label)
            if key == "sp":
                ax.axhline(0, color=MUTED, lw=0.8)
            if target is not None:
                ax.axhline(target, color=MUTED, ls=":", lw=1.1, label=f"target {target:g}")
                ax.legend(loc="upper right")
            e = ex.get(key)
            if e:
                v = e["value"] * (mpa if key in ("vm", "sp") else 1)
                ax.plot([e["x"] * 1e3], [v], "o", color=color, ms=5)
                ax.annotate(f"{v:.2f}: {e['text'].rsplit(', x =', 1)[0]}", xy=(e["x"] * 1e3, v),
                            xytext=(8, -12 if key == "vm" else 4), textcoords="offset points", fontsize=8,
                            color=color)  # the maximum sits at the top of its panel: label below it
        a3.set_xlabel("x [mm]  (0 = throat)")
        self._finish("stress")

    def _draw_components(self, res, sres):
        """Every checked stress on its own: the load components, then the stresses at each check point."""
        x, mpa = res["st"]["x"] * 1e3, 1e-6
        fig, (a1, a2, a3) = self._ax("components", f"Individual stresses ({sres['material'].name}; "
                                                   "positive = tension; fully elastic)", nrows=3)
        a1.plot(x, sres["sb"] * mpa, color=PC["cool"], lw=1.5, label="σ_b,root  pressure bending at rib root")
        a1.plot(x, sres["sb_mid"] * mpa, color=PC["sky"], ls="--", lw=1.3, label="σ_b,mid  pressure bending at mid-span")
        a1.plot(x, sres["st"] * mpa, color=PC["warm"], lw=1.5, label="σ_t  through-wall thermal (magnitude)")
        a1.set_ylabel("tangential loads\n[MPa]")
        a1.legend(loc="upper left")
        a2.plot(x, sres["sl"] * mpa, color=PC["hot"], lw=1.5, label="σ_l  longitudinal thermal (compressive)")
        a2.plot(x, sres["sc"] * 0.9 * mpa, color=PC["green"], ls="--", lw=1.3, label="0.9·S_c  buckling limit (elastic)")
        a2.set_ylabel("longitudinal\n[MPa]")
        a2.set_ylim(bottom=0)
        a2.legend(loc="upper left")
        styles = [(PC["hot"], "-"), (PC["hot"], ":"), (PC["cool"], "-"), (PC["cool"], ":")]
        for p, (face, pos, _, _) in enumerate(ws.POINTS):
            color, ls = styles[p]
            a3.plot(x, sres["vm_p"][:, p] * mpa, color=color, ls=ls, lw=1.4, label=f"σ_vm, {face}, {pos}")
        a3.set_ylabel("von Mises per\ncheck point [MPa]")
        a3.set_xlabel("x [mm]  (0 = throat)")
        a3.legend(loc="upper left", ncol=2)
        self._finish("components")

    def _draw_life(self, res, sres):
        x, ex, si = res["st"]["x"] * 1e3, sres["extremes"], sres["inputs"]
        fig, (a1, a2) = self._ax("life", f"Low-cycle fatigue and longitudinal buckling ({sres['material'].name}; "
                                         "elastic strain range, elastic modulus)",
                                 nrows=2)
        if np.isfinite(sres["nf"]).any():
            a1.semilogy(x, sres["nf"], color=PC["violet"], lw=1.6, label=f"cycles to failure ({sres['life_method']})")
            a1.semilogy(x, sres["nf"] / si.life_factor, color=PC["sky"], ls="--", lw=1.3,
                        label=f"allowable firings (÷{si.life_factor:g})")
            e = ex["nf"]
            a1.plot([e["x"] * 1e3], [e["value"]], "o", color=PC["violet"], ms=5)
            a1.annotate(f"{e['value']:.3g} cycles: {e['text'].rsplit(', x =', 1)[0]}, strain {e['eps'] * 100:.2f} %",
                        xy=(e["x"] * 1e3, e["value"]), xytext=(8, -12), textcoords="offset points", fontsize=8,
                        color=PC["violet"])
            a1.legend(loc="upper left")
        else:
            a1.text(0.5, 0.5, f"{sres['material'].name}: no fatigue data", ha="center", va="center",
                    color=MUTED, transform=a1.transAxes)
        a1.set_ylabel("cycles")
        a2.plot(x, sres["bm"], color=PC["green"], lw=1.6, label="margin 0.9·S_c / σ_l − 1 (elastic S_c)")
        a2.axhline(0, color=BAD, ls=":", lw=1.1, label="zero margin")
        if "bm" in ex:
            e = ex["bm"]
            a2.plot([e["x"] * 1e3], [e["value"]], "o", color=PC["green"], ms=5)
            a2.annotate(f"{e['value']:.2f}", xy=(e["x"] * 1e3, e["value"]), xytext=(8, 4),
                        textcoords="offset points", fontsize=8, color=PC["green"])
        a2.set_ylabel("buckling margin")
        a2.set_xlabel("x [mm]  (0 = throat)")
        a2.legend(loc="lower right")
        self._finish("life")

    def _draw_budget(self, feed):
        b, bar = feed["budget"], 1 / eng.BAR
        fig, ax = self._ax("budget", "Feed pressure budget")
        sides = ["Fuel (through regen)", "Oxidizer"]
        layers = [("chamber pressure", [b["pc"], b["pc"]], PC["bar_pc"]),
                  ("injector ΔP", [b["dp_inj_f"], b["dp_inj_o"]], PC["bar_inj"]),
                  ("regen channel ΔP", [b["dp_regen"], 0.0], PC["bar_regen"])]
        bottom = np.zeros(2)
        top = max(b["p_feed_f"], b["p_feed_o"]) * bar
        for label, vals, color in layers:
            vals = np.array(vals) * bar
            ax.bar(sides, np.where(vals > 0, vals, np.nan), bottom=bottom, color=color, edgecolor=SURFACE,
                   width=0.45, label=label)
            for k, (v, b0) in enumerate(zip(vals, bottom)):
                if v > 0.04 * top:
                    ax.text(k, b0 + v / 2, f"{v:.2f}", ha="center", va="center", fontsize=8, color=FG)
            bottom += vals
        for k, total in enumerate([b["p_feed_f"] * bar, b["p_feed_o"] * bar]):
            ax.text(k, total + 0.01 * top, f"{total:.2f} bar", ha="center", va="bottom", fontweight="bold")
        ax.set_ylabel("pressure [bar]")
        ax.set_ylim(0, top * 1.15)
        ax.grid(axis="x", visible=False)
        ax.legend(loc="upper center", ncol=3)
        self._finish("budget")

    def _draw_sweep(self, design):
        if not design.sweep:
            self._placeholder("sweep", "Enable the O/F sweep on the Engine tab to see this plot.")
            return
        sw = design.sweep
        of = [q["of"] for q in sw]
        fig, ax = self._ax("sweep", "Mixture ratio sweep (delivered performance)")
        ax.plot(of, [q["isp"] for q in sw], "o-", color=PC["hot"], lw=1.6, ms=4, label="Isp (design ambient)")
        ax.plot(of, [q["isp_vac"] for q in sw], "s-", color=PC["violet"], lw=1.6, ms=4, label="Isp (vacuum)")
        ax.axvline(design.perf["of"], color=MUTED, ls="--", lw=1, label=f"design O/F = {design.perf['of']:.2f}")
        ax.set_xlabel("O/F")
        ax.set_ylabel("Isp [s]")
        ax2 = ax.twinx()
        ax2.plot(of, [q["tc"] for q in sw], "^--", color=PC["warm"], lw=1.2, ms=4, label="chamber temperature")
        ax2.set_ylabel("Tc [K]")
        ax2.grid(False)
        ax2.spines["right"].set_visible(True)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="lower center")
        self._finish("sweep")

    # ---------------------------------------------------------------- buttons
    def _apply_preset(self, _=None):
        name = self.vars["preset"].get()
        if name == list(PRESETS)[0]:
            self._reset(keep_preset=True)
            return
        for k, v in PRESETS[name].items():
            self.vars[k].set(v)
        self._set_status(f"Loaded example case: {name}.")

    def _toggle_auto_width(self):
        """Enable the channel inputs that apply: manual, widths from ribs, or optimizer."""
        opt = self.vars["opt_on"].get()
        auto = self.vars["auto_w"].get() and not opt
        self.auto_w_cb.configure(state="disabled" if opt else "normal")
        self.n_ch_entry.configure(state="disabled" if opt else "normal")
        for en in self.width_entries:
            en.configure(state="disabled" if auto or opt else "normal")
        for en in self.height_entries:
            en.configure(state="disabled" if opt else "normal")
        for en in self.rib_entries:
            en.configure(state="normal" if auto else "disabled")
        for en in self.opt_entries:
            en.configure(state="normal" if opt else "disabled")

    def _apply_material(self, _=None):
        self.vars["k_wall"].set(_fmt(MATERIALS[self.vars["material"].get()]))

    def _browse_table(self):
        path = filedialog.askopenfilename(title="Fuel property table",
                                          filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if path:
            self.vars["table"].set(path)

    def _reset(self, keep_preset=False):
        for k, v in self.defaults.items():
            if keep_preset and k == "preset":
                continue
            self.vars[k].set(v)
        self._toggle_auto_width()
        self._set_status("Inputs reset to defaults.")

    def _save_settings(self):
        """Write every input (as typed, 'auto' included) to a JSON file. Outputs are not saved."""
        path = filedialog.asksaveasfilename(title="Save engine settings", defaultextension=".json",
                                            filetypes=[("Engine settings", "*.json"), ("All files", "*.*")])
        if not path:
            return
        data = dict(app="Regen Engine Designer", version=SETTINGS_VERSION,
                    inputs={k: v.get() for k, v in self.vars.items() if k != "preset"})
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except OSError as exc:
            messagebox.showerror("Save settings", f"Could not write {path}:\n{exc}")
            return
        self._set_status(f"Settings saved to {os.path.basename(path)}.")

    def _load_settings(self):
        """Fill the inputs from a settings file. Inputs missing from the file keep their current value."""
        path = filedialog.askopenfilename(title="Load engine settings",
                                          filetypes=[("Engine settings", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as f:
                inputs = json.load(f)["inputs"]
            if not isinstance(inputs, dict):
                raise ValueError("'inputs' is not a table of values")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            messagebox.showerror("Load settings", f"{os.path.basename(path)} is not a settings file:\n{exc}")
            return
        ignored = 0
        for key, val in inputs.items():
            var = self.vars.get(key)
            if var is None or key == "preset":
                ignored += 1
            elif isinstance(var, tk.BooleanVar):
                var.set(bool(val))
            else:
                var.set(str(val))
        missing = sum(1 for k in self.vars if k != "preset" and k not in inputs)
        self._toggle_auto_width()
        msg = f"Settings loaded from {os.path.basename(path)}."
        if ignored or missing:
            msg += f" {ignored} unknown entr{'y' if ignored == 1 else 'ies'} ignored, {missing} input" \
                   f"{'' if missing == 1 else 's'} not in the file kept as they were."
        self._set_status(msg)

    def _export(self):
        if not self.result:
            messagebox.showinfo("Export results", "Run the analysis before exporting.")
            return
        folder = filedialog.askdirectory(title="Export results to folder")
        if not folder:
            return
        design, res, texts = self.result["design"], self.result["res"], self.result["texts"]
        report = texts["engine"] + texts["thermal"] + texts["feed"] + texts.get("stress", "")
        with contextlib.redirect_stdout(io.StringIO()):
            eng.write_csv(design, os.path.join(folder, "engine_contour.csv"))
            if res:
                rt.write_thermal_csv(res, os.path.join(folder, "thermal_profile.csv"))
        with open(os.path.join(folder, "report.txt"), "w", encoding="utf-8") as f:
            f.write(report)
        n = 2 + (1 if res else 0)
        for key, (fig, _) in self.figs.items():
            fig.savefig(os.path.join(folder, f"plot_{key}.png"), dpi=150)
            n += 1
        messagebox.showinfo("Export results", f"Saved {n} files to:\n{folder}")

    def _help(self):
        messagebox.showinfo("Help", (
            "1. Engine: operating point, propellants and mixture ratio (or load an example case).\n"
            "2. Geometry: chamber L*, contraction ratio, cone angles and efficiencies.\n"
            "3. Cooling: fuel property model, inlet temperature and liner wall.\n"
            "4. Channels: channel count and dimensions at the cylinder, throat and exit.\n"
            "5. Injector: stiffness, hole count and discharge coefficient for each side.\n"
            "6. Click Run analysis.\n\n"
            "Summary lists every result and any warnings. The other tabs are plots; use the toolbar "
            "to zoom and pan. Export results saves CSVs, the text report and all plots.\n\n"
            "The fuel is the coolant: it passes through the regen channels and then the injector. "
            "Its inlet pressure is solved so it reaches the injector at Pc × (1 + stiffness).\n\n"
            "Hover over any input for a description. 'auto' lets the code choose the value."))

    def _close(self):
        plt.close("all")
        self.quit()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
