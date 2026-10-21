"""
POOP INDUSTRIES(TM) Regen Engine Designer 9000 Pro Max: Deep Fried Edition.

A meme skin over rocket_gui.App. Same inputs, same physics, same results; only the colours, fonts,
labels and vibes change. The real GUI (rocket_gui.py) is untouched and keeps its own theme setting.

    py -3.14 rocket_gui_meme.py
"""

from __future__ import annotations

import os
import random
import tkinter as tk
from tkinter import ttk

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

import rocket_gui as rg

# --------------------------------------------------------------------------------------
# Fonts: Comic Sans for everything, Impact for anything important (as is tradition)
# --------------------------------------------------------------------------------------
rg.UI = "Comic Sans MS"
rg.UI_DISPLAY = "Impact"
rg.F = (rg.UI, 10)
rg.FB = (rg.UI, 10, "bold")
rg.MONO = ("Lucida Console", 9)        # reports stay monospaced so the tables still line up

# --------------------------------------------------------------------------------------
# Palettes. "light" = Normie Mode (fresh poop), "dark" = Deep Fried Mode (well-done poop).
# Every value unique within a palette (rocket_gui repaints by role).
# --------------------------------------------------------------------------------------
rg.THEMES["light"] = dict(bg="#f2dfb3", surface="#fff3d1", border="#b98550", fg="#3a1d0b", muted="#7b4a24",
                          accent="#ff00cc", ok="#00a400", warn="#e07000", bad="#e00000")
rg.THEMES["dark"] = dict(bg="#24130a", surface="#3a1f0d", border="#7a4518", fg="#ffe14d", muted="#e0a061",
                         accent="#39ff14", ok="#7cfc00", warn="#ff9f1a", bad="#ff2e63")
rg.PLOT_COLORS["light"] = dict(grid="#ead7ad", edge="#8b5a2b", hot="#ff0000", warm="#ff8c00", cool="#0066ff",
                               aw="#a0522d", violet="#9400d3", green="#00a400", sky="#00bfff", gas="#ffd27f",
                               liner="#6b3a1e", channel="#4fa3ff", bar_pc="#c8a26b", bar_inj="#ff69b4",
                               bar_regen="#8b4513")
rg.PLOT_COLORS["dark"] = dict(grid="#4a2a12", edge="#c08040", hot="#ff3b3b", warm="#ffb000", cool="#00e5ff",
                              aw="#c8875a", violet="#ff66ff", green="#39ff14", sky="#7fdbff", gas="#5a3415",
                              liner="#a0602a", channel="#1e90ff", bar_pc="#6b4423", bar_inj="#ff1493",
                              bar_regen="#d2691e")

# Own preferences file, so the meme theme choice never touches the real GUI's setting
rg.PREFS_PATH = os.path.join(os.path.expanduser("~"), ".poop_industries.json")

# Hand-drawn xkcd plots, the most 2012 thing matplotlib can do
plt.xkcd(scale=0.8, length=120, randomness=2)
_real_use_palette = rg._use_palette


def _use_palette(name: str) -> None:
    _real_use_palette(name)
    # xkcd draws a halo round every line and letter; make it the plot background, not white
    plt.rcParams["path.effects"] = [pe.withStroke(linewidth=3, foreground=rg.SURFACE)]


rg._use_palette = _use_palette
_use_palette(rg._load_theme())

# --------------------------------------------------------------------------------------
# Words
# --------------------------------------------------------------------------------------
BRAND = "POOP INDUSTRIES™"
SLOGANS = [
    "we put the 'turd' in 'turbopump'",
    "proudly #2 in rocketry",
    "engineered in a toilet, tested in a toilet",
    "now with 30% more brown",
    "it's not a leak, it's film cooling",
    "the only engine that runs on regret",
    "certified fresh by nobody",
    "where every design is a number two",
]
DOGE = ["such regen", "very cool", "much heat", "so throat", "wow", "many channel", "very nozzle",
        "such flux", "much pressure", "so stress", "wow such Isp"]
STATUS_SUFFIX = ["💩", "(no cap)", "#YOLO", "#swag", "lol", "xD", "( ͡° ͜ʖ ͡°)", "¯\\_(ツ)_/¯", "🔥🔥🔥",
                 "rawr XD", "le epic", "*dabs*", "much wow", "420 blaze it", "#blessed", "Harlem Shake!!"]
TICKER = ("💩 POOP INDUSTRIES™ 💩  ·  ONE DOES NOT SIMPLY COOL A THROAT  ·  Y U NO CONVERGE  ·  "
          "IT'S OVER 9000 (KELVIN)  ·  SHUT UP AND TAKE MY THRUST  ·  I CAN HAZ ISP?  ·  "
          "ERMAHGERD, RERKERT  ·  DO A BARREL ROLL  ·  ALL YOUR BASE ARE BELONG TO REGEN  ·  "
          "NOT SURE IF FILM COOLING OR JUST LEAKING  ·  THE CAKE IS A LIE BUT THE HEAT FLUX IS NOT  ·  "
          "CALL ME MAYBE (WHEN THE WALL MELTS)  ·  ")
TABS = {
    "Engine": "Engine 🚀", "Geometry": "Shape ™", "Cooling": "Cool Zone 😎", "Channels": "Poop Chutes 💩",
    "Injector": "Squirt Holes 💦", "Summary": "The Tea ☕", "Structures": "Will It Break? 😬",
    "Structural plots": "Stress Graphs 📉", "Temperatures": "Hot Stuff 🔥",
    "Heat flux": "Flux Capacitor ⚡", "Coolant": "Poop Juice 🟤", "Pressure budget": "Pressure 😤",
    "O/F sweep": "O/F Sweep 🧹", "Envelope & yield SF": "Envelope & Yield SF (yikes)",
    "Stress components": "Stress Bits", "Fatigue & buckling": "Fatigue & Buckling (big oof)",
}
CARDS = {
    "Isp (design)": "Isp (design) 🚀", "Isp (vacuum)": "Isp (space lol)", "Thrust": "Yeet Force",
    "Mass flow": "Poop Flow Rate", "c*": "c* (c-star, not c-poop)", "Throat diameter": "Throat Goat",
    "Exit diameter": "Butthole Diameter", "Max hot-wall T": "Hottest Wall 🔥", "Regen ΔP": "Regen ΔP 😤",
    "Feed pressure": "Squeeze Pressure", "Orifice diameter": "Squirt Hole Size",
    "Liquid film length": "Skid Mark Length", "Max von Mises stress": "von Mises (big stress)",
    "Max principal stress": "Principal Stress (the boss)", "Min yield SF": "Yield SF (yikes?)",
    "Fatigue life": "Fatigue Life (tired)", "Buckling margin": "Buckle Up Margin",
}
ROLE_EMOJI = {"ok": " 👌", "warn": " 😬", "bad": " 💀"}


class MemeApp(rg.App):
    def __init__(self):
        super().__init__()
        self.title(f"💩 {BRAND} Regen Engine Designer 9000 Pro Max: Deep Fried Edition 💩")
        self.go_btn.configure(text="💩 SEND IT 💩")
        self._rename_tabs(self)
        self._ticker()
        self._say(f"Welcome to {BRAND}.\n\nSet the inputs on the left (or load an example, they're good I promise) "
                  "and smash that SEND IT button like it's 2012.\n\nDon't forget to like and subscribe.", tag="head")
        self._set_status("Hover over anything for a description. Or don't. I'm a GUI, not a cop.")

    # ---------------------------------------------------------------- header + ticker
    def _header(self):
        top = self._reg(tk.Frame(self), bg="bg")
        top.pack(fill="x")
        txt = self._reg(tk.Frame(top), bg="bg")
        txt.pack(side="left", padx=18, pady=(10, 4))
        row = self._reg(tk.Frame(txt), bg="bg")
        row.pack(anchor="w")
        self._reg(tk.Label(row, text=f"💩 {BRAND}", font=("Impact", 26)), bg="bg", fg="accent").pack(side="left")
        self._reg(tk.Label(row, text="  presents", font=("Comic Sans MS", 11, "italic")),
                  bg="bg", fg="muted").pack(side="left", anchor="s", pady=(0, 6))
        self._reg(tk.Label(txt, text="REGEN ENGINE DESIGNER 9000 PRO MAX (DEEP FRIED EDITION)",
                           font=("Impact", 15)), bg="bg", fg="fg").pack(anchor="w")
        self.slogan = self._reg(tk.Label(txt, text=f"“{random.choice(SLOGANS)}”", font=("Comic Sans MS", 10)),
                                bg="bg", fg="muted")
        self.slogan.pack(anchor="w")
        ttk.Button(top, text="Halp 🆘", command=self._help).pack(side="right", padx=(6, 18))
        ttk.Checkbutton(top, text="Deep fried mode 🔥", variable=self.dark_var, command=self._toggle_theme,
                        style="Switch.TCheckbutton" if rg.sv_ttk else "TCheckbutton").pack(side="right", padx=6)
        self.verdict = self._reg(tk.Label(top, text="", font=("Impact", 18)), bg="bg", fg="accent")
        self.verdict.pack(side="right", padx=20)
        bar = self._reg(tk.Frame(self, height=self.px(26)), bg="border")
        bar.pack(fill="x")
        bar.pack_propagate(False)
        self.ticker = self._reg(tk.Label(bar, text=TICKER * 3, font=("Impact", 11), anchor="w"), bg="border",
                                fg="fg")
        self.ticker.place(x=0, rely=0.5, anchor="w")
        self._tick_x = 0

    def _ticker(self):
        """Scroll the marquee like a 2009 GeoCities page."""
        try:
            width = self.ticker.winfo_reqwidth() // 3
            self._tick_x = (self._tick_x - 2) % max(width, 1)
            self.ticker.place_configure(x=-self._tick_x)
            self.after(30, self._ticker)
        except tk.TclError:  # window closed
            pass

    def _style(self):
        super()._style()
        s = ttk.Style(self)
        s.configure("Accent.TButton", font=("Impact", 16))
        s.configure("Run.TButton", font=("Impact", 16), padding=(10, 8))
        s.configure("TNotebook.Tab", font=("Comic Sans MS", 10, "bold"))

    # ---------------------------------------------------------------- words everywhere
    def _rename_tabs(self, widget):
        for child in widget.winfo_children():
            if isinstance(child, ttk.Notebook):
                for tab in child.tabs():
                    old = child.tab(tab, "text")
                    child.tab(tab, text=TABS.get(old, old))
            self._rename_tabs(child)

    def _section(self, parent, text, **kw):
        return super()._section(parent, f"{text} {random.choice(['💩', '🚀', '🔥', '😎', '👌', '💯', '🤔'])}", **kw)

    def _cards(self, parent, names):
        """Same tiles as the real GUI (keyed by the real names), with better (worse) labels."""
        cards = ttk.Frame(parent)
        cards.pack(fill="x", pady=(0, 6))
        for i, name in enumerate(names):
            c = self._reg(tk.Frame(cards, highlightthickness=2, padx=14, pady=10),
                          bg="surface", highlightbackground="border")
            c.grid(row=i // 6, column=i % 6, sticky="nsew", padx=3, pady=3)
            cards.columnconfigure(i % 6, weight=1, uniform="cards")
            self._reg(tk.Label(c, text=CARDS.get(name, name), font=("Comic Sans MS", 9, "bold"), anchor="w"),
                      bg="surface", fg="muted").pack(fill="x")
            val = self._reg(tk.Label(c, text="—", font=("Impact", 17), anchor="w"), bg="surface", fg="fg")
            val_roles = self._themed[-1][1]
            val.pack(fill="x")
            sub = self._reg(tk.Label(c, text="", font=("Comic Sans MS", 8), anchor="w"), bg="surface", fg="muted")
            sub.pack(fill="x")
            self.cards[name] = (val, sub, val_roles)

    def _card(self, name, value, sub, role="fg"):
        super()._card(name, value + ROLE_EMOJI.get(role, ""), sub, role)

    def _set_status(self, text, role="muted"):
        super()._set_status(f"{text}  {random.choice(STATUS_SUFFIX)}", role)

    def _placeholder(self, key, msg="Smash SEND IT to see results (no cap) 💩"):
        super()._placeholder(key, msg)

    def _ax(self, key, title, nrows=1, **kw):
        fig, axes = super()._ax(key, f"{title}   ({random.choice(DOGE)})", nrows, **kw)
        fig.text(0.995, 0.005, f"© {BRAND} — all rights reserved, no refunds", ha="right", va="bottom",
                 fontsize=8, color=rg.MUTED, alpha=0.8)
        return fig, axes

    # ---------------------------------------------------------------- run
    def run(self):
        super().run()
        if self.busy:
            self.go_btn.configure(text=random.choice(["LOADING… 💩", "DOING A BARREL ROLL…", "COMPUTING THE POOPS…",
                                                      "HOLD MY BEER…", "PLEASE WAIT (OR DON'T)…"]))
            self.verdict.configure(text="")

    def _poll(self):
        super()._poll()
        if self.busy:
            return
        self.go_btn.configure(text="💩 SEND IT 💩")
        self.slogan.configure(text=f"“{random.choice(SLOGANS)}”")
        if self.result and self.result.get("res"):
            t_max = self._tmax(self.result["res"])
            limit = self.result["th"].t_wall_limit_K
            self.verdict.configure(text="THIS IS FINE 🔥🐶" if t_max > limit else "SUCH COOL. WOW. 😎")
        elif self.result:
            self.verdict.configure(text="NO COOLING?? BOLD.")


if __name__ == "__main__":
    MemeApp().mainloop()
