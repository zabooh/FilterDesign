"""FilterDesign GUI - interactive IIR filter design with C code generation.

Run:  python app.py  [--host 127.0.0.1] [--port 8090] [--no-browser]

The design is computed by the fdesign command line tool (../fdesign), the
analysis and the code generation are done in Python (fdcore.py, codegen.py).
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from nicegui import ui
from scipy import signal

import codegen
import fdcore

TYPE_LABELS = {"lowpass": "Lowpass", "highpass": "Highpass", "bandpass": "Bandpass", "bandstop": "Bandstop"}
CHAR_LABELS = {"butterworth": "Butterworth", "chebyshev": "Chebyshev (type I)", "elliptic": "Elliptic (Cauer)"}
ARITH_LABELS = {"float": "float (32 bit)", "double": "double (64 bit)",
                "fixed32": "fixed point 32 bit", "fixed16": "fixed point 16 bit"}

# edge widgets in ascending frequency order: (label, list, index)
EDGE_LAYOUT = {
    "lowpass": [("Passband edge", "fpass", 0), ("Stopband edge", "fstop", 0)],
    "highpass": [("Stopband edge", "fstop", 0), ("Passband edge", "fpass", 0)],
    "bandpass": [("Lower stopband edge", "fstop", 0), ("Lower passband edge", "fpass", 0),
                 ("Upper passband edge", "fpass", 1), ("Upper stopband edge", "fstop", 1)],
    "bandstop": [("Lower passband edge", "fpass", 0), ("Lower stopband edge", "fstop", 0),
                 ("Upper stopband edge", "fstop", 1), ("Upper passband edge", "fpass", 1)],
}
# default edges as fraction of fs, ascending
EDGE_DEFAULTS = {"lowpass": [0.05, 0.1], "highpass": [0.15, 0.225],
                 "bandpass": [0.05, 0.1, 0.15, 0.225], "bandstop": [0.05, 0.1, 0.15, 0.225]}

N_FREQ = 2048
N_TIME_MAX = 4096
PALETTES = {
    "light": {"design": "#2563eb", "impl": "#f59e0b", "zero": "#2563eb", "pole": "#dc2626",
              "forbidden": "rgba(220, 38, 38, 0.10)", "grid": "#e5e7eb", "fg": "#374151",
              "axis": "#6b7280", "circle": "#9ca3af"},
    "dark": {"design": "#60a5fa", "impl": "#fbbf24", "zero": "#60a5fa", "pole": "#f87171",
             "forbidden": "rgba(248, 113, 113, 0.14)", "grid": "#334155", "fg": "#cbd5e1",
             "axis": "#64748b", "circle": "#64748b"},
}
START_DARK = True


def r6(x: float) -> float | None:
    if x is None or not math.isfinite(x):
        return None
    return float(f"{x:.6g}")


def xy(x: np.ndarray, y: np.ndarray) -> list:
    return [[r6(a), r6(b)] for a, b in zip(x, y)]


def base_chart(xname: str, yname: str, xtype: str = "value") -> dict:
    return {
        "animation": False,
        "grid": {"left": 70, "right": 30, "top": 40, "bottom": 80},
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "cross"}},
        "legend": {"top": 5},
        "toolbox": {"right": 10, "feature": {"dataZoom": {}, "restore": {}, "saveAsImage": {}}},
        "xAxis": {"type": xtype, "name": xname, "nameLocation": "middle", "nameGap": 28},
        "yAxis": {"type": "value", "name": yname, "nameLocation": "middle", "nameGap": 50},
        "dataZoom": [{"type": "inside", "xAxisIndex": 0, "filterMode": "none"},
                     {"type": "slider", "xAxisIndex": 0, "filterMode": "none", "height": 18, "bottom": 10}],
        "series": [],
    }


def apply_theme(opt: dict, pal: dict) -> dict:
    """Colors text, axes, grid lines and controls of an ECharts option dict."""
    fg = {"color": pal["fg"]}
    opt["backgroundColor"] = "transparent"
    opt["textStyle"] = fg
    opt.setdefault("legend", {})["textStyle"] = fg
    for ax in ("xAxis", "yAxis"):
        a = opt.setdefault(ax, {})
        a["axisLine"] = {"lineStyle": {"color": pal["axis"]}}
        a["axisLabel"] = fg
        a["nameTextStyle"] = fg
        a["splitLine"] = {"lineStyle": {"color": pal["grid"]}}
    if "toolbox" in opt:
        opt["toolbox"]["iconStyle"] = {"borderColor": pal["fg"]}
    for dz in opt.get("dataZoom", []):
        if dz.get("type") == "slider":
            dz.update({"textStyle": fg, "borderColor": pal["axis"], "fillerColor": "rgba(96,165,250,0.15)",
                       "dataBackground": {"lineStyle": {"color": pal["axis"]}, "areaStyle": {"color": pal["grid"]}}})
    if "tooltip" in opt:
        opt["tooltip"].update({"backgroundColor": "#1e293b" if pal is PALETTES["dark"] else "#ffffff",
                               "borderColor": pal["axis"], "textStyle": fg})
    return opt


def line(name: str, data: list, color: str, dashed: bool = False, width: float = 2) -> dict:
    return {"name": name, "type": "line", "data": data, "showSymbol": False, "connectNulls": False,
            "lineStyle": {"color": color, "width": width, "type": "dashed" if dashed else "solid"},
            "itemStyle": {"color": color}}


TAB_NAMES = ["Magnitude", "Phase", "Group delay", "Pole / zero", "Impulse / step", "Sections", "C code"]


@ui.page("/")
def index(tab: str = "magnitude") -> None:
    """Main page; ?tab=phase etc. opens a specific tab."""
    dark = ui.dark_mode(value=START_DARK)

    def P() -> dict:
        return PALETTES["dark" if dark.value else "light"]

    state = {
        "spec": fdcore.Spec(),
        "edges": {t: None for t in EDGE_LAYOUT},
        "design": None,
        "timer": None,
    }
    # initial bandpass edges like the firmware example
    state["edges"]["bandpass"] = [400.0, 800.0, 1200.0, 1800.0]

    ui.add_head_html("<style>.q-field--dense .q-field__control{height:36px}</style>")

    # ------------------------------------------------------------------ header
    with ui.header().classes("items-center bg-slate-800 text-white px-4 py-2"):
        ui.label("FilterDesign").classes("text-xl font-semibold")
        ui.label("IIR filter design · Butterworth · Chebyshev · Elliptic").classes("text-sm opacity-70 ml-3")
        ui.space()
        theme_btn = ui.button(icon="light_mode").props("flat round color=white").tooltip("Light / dark mode")
        status_chip = ui.label("").classes("text-sm px-3 py-1 rounded")

    with ui.row().classes("w-full no-wrap items-start gap-4 p-4"):
        # ============================================================== left panel
        with ui.column().classes("w-80 shrink-0 gap-3"):
            with ui.card().classes("w-full"):
                ui.label("Specification").classes("text-base font-semibold")
                type_sel = ui.select(TYPE_LABELS, value="bandpass", label="Filter type").props("dense").classes("w-full")
                char_sel = ui.select(CHAR_LABELS, value="elliptic", label="Characteristic").props("dense").classes("w-full")
                fs_in = ui.number("Sampling rate", value=8000, min=1, suffix="Hz", format="%.6g").props("dense").classes("w-full")
                edge_in = [ui.number(f"Edge {i + 1}", value=0, min=0, suffix="Hz", format="%.6g").props("dense").classes("w-full")
                           for i in range(4)]
                with ui.row().classes("w-full no-wrap gap-2"):
                    ap_in = ui.number("Passband ripple", value=1.0, min=0.001, step=0.1, suffix="dB",
                                      format="%.4g").props("dense").classes("w-1/2")
                    as_in = ui.number("Stopband atten.", value=40.0, min=1, step=1, suffix="dB",
                                      format="%.4g").props("dense").classes("w-1/2")

            with ui.card().classes("w-full"):
                ui.label("Implementation").classes("text-base font-semibold")
                arith_sel = ui.select(ARITH_LABELS, value="fixed32", label="Arithmetic").props("dense").classes("w-full")
                scale_chk = ui.checkbox("Section scaling (L∞, 0 dB per node)", value=True)
                with ui.row().classes("w-full no-wrap items-center gap-2"):
                    auto_q = ui.checkbox("Auto Q", value=True)
                    frac_in = ui.number("Fractional bits", value=30, min=0, max=31, step=1,
                                        format="%d").props("dense").classes("grow")
                name_in = ui.input("C name", value="iir_filter").props("dense").classes("w-full")

            error_box = ui.label("").classes("w-full text-sm text-red-500 bg-red-500/10 p-2 rounded hidden")

        # ============================================================== right panel
        with ui.column().classes("grow min-w-0"):
            with ui.card().classes("w-full py-2"):
                summary = ui.grid(columns=3).classes("w-full gap-x-10 gap-y-0 text-sm")
                warn_box = ui.column().classes("w-full gap-0 text-sm")
            with ui.tabs().classes("w-full") as tabs:
                t_mag = ui.tab("Magnitude")
                t_phase = ui.tab("Phase")
                t_gd = ui.tab("Group delay")
                t_pz = ui.tab("Pole / zero")
                t_time = ui.tab("Impulse / step")
                t_sec = ui.tab("Sections")
                t_code = ui.tab("C code")
            tab_widgets = [t_mag, t_phase, t_gd, t_pz, t_time, t_sec, t_code]
            wanted = next((w for w, n in zip(tab_widgets, TAB_NAMES)
                           if n.lower().replace(" ", "").replace("/", "").startswith(tab.lower())), t_mag)
            with ui.tab_panels(tabs, value=wanted).classes("w-full"):
                with ui.tab_panel(t_mag):
                    with ui.row().classes("items-center gap-6"):
                        mag_scale = ui.toggle({"db": "dB", "lin": "linear"}, value="db")
                        mag_xlog = ui.toggle({"lin": "linear f", "log": "log f"}, value="lin")
                        mag_mask = ui.checkbox("Tolerance scheme", value=True)
                        mag_impl = ui.checkbox("Implemented (quantized)", value=True)
                    mag_chart = ui.echart(base_chart("Frequency (Hz)", "Magnitude (dB)")).classes("w-full h-[520px]")
                with ui.tab_panel(t_phase):
                    phase_wrap = ui.toggle({"unwrap": "unwrapped", "wrap": "wrapped (±180°)"}, value="unwrap")
                    phase_chart = ui.echart(base_chart("Frequency (Hz)", "Phase (°)")).classes("w-full h-[520px]")
                with ui.tab_panel(t_gd):
                    gd_unit = ui.toggle({"samples": "samples", "ms": "ms"}, value="samples")
                    gd_chart = ui.echart(base_chart("Frequency (Hz)", "Group delay (samples)")).classes("w-full h-[520px]")
                with ui.tab_panel(t_pz):
                    with ui.row().classes("w-full items-start gap-6"):
                        pz_chart = ui.echart({"series": []}).classes("w-[560px] h-[560px]")
                        pz_table = ui.column().classes("text-sm")
                with ui.tab_panel(t_time):
                    time_kind = ui.toggle({"impulse": "impulse response", "step": "step response"}, value="impulse")
                    time_chart = ui.echart(base_chart("Sample n", "Amplitude")).classes("w-full h-[520px]")
                with ui.tab_panel(t_sec):
                    sec_info = ui.label("").classes("text-sm opacity-70")
                    sec_table = ui.table(columns=[], rows=[], row_key="k").classes("w-full").props("dense flat")
                with ui.tab_panel(t_code):
                    with ui.row().classes("items-center gap-3"):
                        dl_h = ui.button("Download .h", icon="download")
                        dl_c = ui.button("Download .c", icon="download")
                        code_note = ui.label("").classes("text-sm opacity-70")
                    with ui.row().classes("w-full no-wrap gap-4 items-start"):
                        with ui.column().classes("w-1/2 min-w-0"):
                            h_title = ui.label("").classes("font-mono text-sm")
                            h_code = ui.code("", language="c").classes("w-full text-xs")
                        with ui.column().classes("w-1/2 min-w-0"):
                            c_title = ui.label("").classes("font-mono text-sm")
                            c_code = ui.code("", language="c").classes("w-full text-xs")

    code_files = {"h": ("", ""), "c": ("", "")}
    dl_h.on_click(lambda: ui.download.content(code_files["h"][1], code_files["h"][0]))
    dl_c.on_click(lambda: ui.download.content(code_files["c"][1], code_files["c"][0]))

    # ------------------------------------------------------------------ helpers
    def load_edges() -> None:
        t = type_sel.value
        if state["edges"][t] is None:
            state["edges"][t] = [round(f * (fs_in.value or 8000), 6) for f in EDGE_DEFAULTS[t]]
        layout = EDGE_LAYOUT[t]
        for i, w in enumerate(edge_in):
            if i < len(layout):
                w.props(f'label="{layout[i][0]}"')
                w.set_value(state["edges"][t][i])
                w.set_visibility(True)
            else:
                w.set_visibility(False)

    def read_spec() -> fdcore.Spec:
        t = type_sel.value
        layout = EDGE_LAYOUT[t]
        vals = [float(edge_in[i].value or 0) for i in range(len(layout))]
        state["edges"][t] = vals
        fpass, fstop = [0.0, 0.0], [0.0, 0.0]
        for (label, lst, idx), v in zip(layout, vals):
            (fpass if lst == "fpass" else fstop)[idx] = v
        return fdcore.Spec(type=t, characteristic=char_sel.value, fs=float(fs_in.value or 0),
                           fpass=fpass, fstop=fstop, ap=float(ap_in.value or 0), as_=float(as_in.value or 0))

    def implementation(d: fdcore.Design):
        sos = fdcore.scale_sections(d.sos, d.spec.fs) if scale_chk.value else d.sos.copy()
        arith = arith_sel.value
        fixed = None
        if arith.startswith("fixed"):
            word = int(arith[5:])
            frac = fdcore.auto_frac_bits(sos, word)
            if auto_q.value:
                frac_in.set_value(frac)
            else:
                frac = int(max(0, min(word - 1, frac_in.value or 0)))
            fixed = fdcore.quantize(sos, word, frac)
            impl_sos = fixed.as_float_sos()
        elif arith == "float":
            impl_sos = sos.astype(np.float32).astype(float)
        else:
            impl_sos = sos
        return sos, fixed, impl_sos

    def schedule(*_):
        if state["timer"] is not None:
            state["timer"].cancel()
        state["timer"] = ui.timer(0.3, update, once=True)

    # ------------------------------------------------------------------ main update
    def update() -> None:
        state["timer"] = None
        try:
            spec = read_spec()
            d = fdcore.design(spec)
        except fdcore.DesignError as e:
            error_box.text = str(e)
            error_box.classes(remove="hidden")
            status_chip.text = "design failed"
            status_chip.classes(replace="text-sm px-3 py-1 rounded bg-red-600")
            return
        error_box.classes(add="hidden")
        state["design"] = d
        fs = spec.fs
        sos, fixed, impl_sos = implementation(d)

        # ---------------- frequency domain
        f = fdcore.freq_axis(fs, N_FREQ, log=mag_xlog.value == "log")
        h = fdcore.response(d.sos, f, fs)
        hi = fdcore.response(impl_sos, f, fs)
        render_magnitude(spec, f, h, hi)
        flin = fdcore.freq_axis(fs, N_FREQ)
        hl = fdcore.response(d.sos, flin, fs)
        render_phase(flin, hl)
        gd = fdcore.group_delay(d.sos, flin, fs)
        render_group_delay(spec, flin, gd)

        # ---------------- poles / zeros
        z, p = fdcore.zeros_poles(d.sos)
        render_pz(z, p, fs)

        # ---------------- time domain
        settle = fdcore.settle_length(d.sos)
        n_time = int(min(max(settle * 1.2, 64), N_TIME_MAX))
        sat = None
        fixed_imp = None
        fixed_time = None
        if fixed is not None:
            fixed_imp, sat = fdcore.fixed_impulse_response(fixed, n_time)
            fixed_time = fdcore.fixed_step_response(fixed, n_time) if time_kind.value == "step" else fixed_imp
        render_time(d, impl_sos, fixed_time, n_time)

        # ---------------- sections, implementation check
        render_sections(d, sos, fixed, fs)
        att_impl_pass = [-20 * math.log10(abs(fdcore.response(impl_sos, [x], fs)[0]) + 1e-300)
                         for x in spec.fpass[:len(d.att_pass)]]
        att_impl_stop = [-20 * math.log10(abs(fdcore.response(impl_sos, [x], fs)[0]) + 1e-300)
                         for x in spec.fstop[:len(d.att_stop)]]
        impl_ok = all(a <= spec.ap + 0.01 for a in att_impl_pass) and all(a >= spec.as_ - 0.01 for a in att_impl_stop)

        # ---------------- summary
        summary.clear()
        pole_r = float(np.max(np.abs(p))) if len(p) else 0.0
        fc = (math.sqrt(spec.fpass[0] * spec.fpass[1]) if spec.type == "bandpass"
              else spec.fpass[0] / 2 if spec.type == "lowpass" else 0.0)
        gd_fc = float(fdcore.group_delay(d.sos, np.array([fc]), fs)[0]) if spec.type in ("lowpass", "bandpass") else None
        fmt = lambda v: ", ".join(f"{x:.2f}" for x in v)
        rows = [
            ("Order", f"{d.digital_order}  ({len(d.sos)} sections)"),
            ("Passband atten.", f"{fmt(d.att_pass)} dB  (≤ {spec.ap:g})"),
            ("Stopband atten.", f"{fmt(d.att_stop)} dB  (≥ {spec.as_:g})"),
            ("Max. pole radius", f"{pole_r:.6f}"),
            ("Settling (−80 dB)", f"{settle} samples = {settle / fs * 1e3:.3g} ms"),
        ]
        if gd_fc is not None and math.isfinite(gd_fc):
            rows.append(("Group delay @ " + f"{fc:.4g} Hz", f"{gd_fc:.2f} samples = {gd_fc / fs * 1e3:.3g} ms"))
        if fixed is not None:
            rows.append(("Coefficient format", f"Q{fixed.word - 1 - fixed.frac}.{fixed.frac} (int{fixed.word}_t)"
                         + ("  CLIPPED!" if fixed.clipped else "")))
            rows.append(("Impulse test (½ FS)", f"{sat} saturation(s)"))
        rows.append(("Impl. passband", f"{fmt(att_impl_pass)} dB"))
        rows.append(("Impl. stopband", f"{fmt(att_impl_stop)} dB"))
        with summary:
            for k, v in rows:
                with ui.row().classes("w-full no-wrap justify-between gap-2"):
                    ui.label(k).classes("opacity-70 whitespace-nowrap")
                    ui.label(v).classes("text-right font-mono whitespace-nowrap")
        warn_box.clear()
        with warn_box:
            if fixed is not None and (fixed.clipped or sat):
                ui.label("Fixed point overflow: enable section scaling or reduce the fractional bits."
                         ).classes("text-red-500 mt-1")
            if not impl_ok:
                ui.label("The quantized implementation violates the specification."
                         ).classes("text-amber-500 mt-1")

        ok = d.spec_met and impl_ok and not (fixed is not None and (fixed.clipped or sat))
        status_chip.text = (f"{CHAR_LABELS[spec.characteristic]} {TYPE_LABELS[spec.type].lower()}, "
                            f"order {d.digital_order} · " + ("specification met" if ok else "check result"))
        status_chip.classes(replace="text-sm px-3 py-1 rounded " + ("bg-green-600" if ok else "bg-amber-600"))

        # ---------------- code
        files = codegen.generate(d, name_in.value or "iir_filter", arith_sel.value, sos, fixed)
        code_files["h"] = (files[0], files[1])
        code_files["c"] = (files[2], files[3])
        h_title.text, c_title.text = files[0], files[2]
        h_code.content, c_code.content = files[1], files[3]
        code_note.text = (f"{ARITH_LABELS[arith_sel.value]}, {len(sos)} sections"
                          + (", section scaling" if scale_chk.value else ""))

    # ------------------------------------------------------------------ renderers
    def show(chart, opt: dict) -> None:
        chart.options.clear()
        chart.options.update(apply_theme(opt, P()))
        chart.update()

    def render_magnitude(spec: fdcore.Spec, f, h, hi) -> None:
        lin = mag_scale.value == "lin"
        opt = base_chart("Frequency (Hz)", "Magnitude" if lin else "Magnitude (dB)",
                         "log" if mag_xlog.value == "log" else "value")
        opt["xAxis"].update({"min": r6(f[0]), "max": r6(f[-1])})
        if lin:
            y, yi = np.abs(h), np.abs(hi)
            opt["yAxis"].update({"min": 0, "max": 1.1})
            pass_lim, stop_lim, top, bottom = 10 ** (-spec.ap / 20), 10 ** (-spec.as_ / 20), 1.1, 0
        else:
            y, yi = fdcore.magnitude_db(h), fdcore.magnitude_db(hi)
            bottom = -10 * math.ceil((spec.as_ + 40) / 10)
            opt["yAxis"].update({"min": bottom, "max": 5})
            pass_lim, stop_lim, top = -spec.ap, -spec.as_, 5
        series = line("Design", xy(f, y), P()["design"])
        if mag_mask.value:
            nyq = spec.fs / 2
            f0 = r6(f[0])
            areas = []
            pa = {"lowpass": [(f0, spec.fpass[0])], "highpass": [(spec.fpass[0], nyq)],
                  "bandpass": [(spec.fpass[0], spec.fpass[1])],
                  "bandstop": [(f0, spec.fpass[0]), (spec.fpass[1], nyq)]}[spec.type]
            sa = {"lowpass": [(spec.fstop[0], nyq)], "highpass": [(f0, spec.fstop[0])],
                  "bandpass": [(f0, spec.fstop[0]), (spec.fstop[1], nyq)],
                  "bandstop": [(spec.fstop[0], spec.fstop[1])]}[spec.type]
            for a, b in pa:
                areas.append([{"xAxis": a, "yAxis": bottom}, {"xAxis": b, "yAxis": pass_lim}])
            for a, b in sa:
                areas.append([{"xAxis": a, "yAxis": stop_lim}, {"xAxis": b, "yAxis": top}])
            series["markArea"] = {"silent": True, "itemStyle": {"color": P()["forbidden"]}, "data": areas}
        opt["series"].append(series)
        if mag_impl.value:
            opt["series"].append(line(f"Implemented ({ARITH_LABELS[arith_sel.value]})", xy(f, yi),
                                      P()["impl"], dashed=True, width=1.5))
        show(mag_chart, opt)

    def render_phase(f, h) -> None:
        ph = fdcore.phase_deg(h) if phase_wrap.value == "unwrap" else np.degrees(np.angle(h))
        opt = base_chart("Frequency (Hz)", "Phase (°)")
        opt["xAxis"].update({"min": 0, "max": r6(f[-1])})
        opt["series"].append(line("Phase", xy(f, ph), P()["design"]))
        show(phase_chart, opt)

    def render_group_delay(spec, f, gd) -> None:
        ms = gd_unit.value == "ms"
        y = gd / spec.fs * 1e3 if ms else gd
        opt = base_chart("Frequency (Hz)", "Group delay (ms)" if ms else "Group delay (samples)")
        opt["xAxis"].update({"min": 0, "max": r6(f[-1])})
        finite = y[np.isfinite(y)]
        if len(finite):
            hi = float(np.percentile(finite, 99.5))
            lo = float(np.min(finite))
            opt["yAxis"].update({"min": min(0, math.floor(lo)), "max": math.ceil(hi * 1.1) if hi > 0 else 1})
        opt["series"].append(line("Group delay", xy(f, y), P()["design"]))
        show(gd_chart, opt)

    def cluster(points: np.ndarray) -> list:
        groups: list[list] = []
        for c in points:
            for g in groups:
                if abs(g[0] - c) < 1e-5:
                    g[1] += 1
                    break
            else:
                groups.append([c, 1])
        return groups

    def render_pz(z, p, fs) -> None:
        rmax = max(1.15, float(np.max(np.abs(np.concatenate([z, p])))) * 1.1 if len(z) + len(p) else 1.15)
        th = np.linspace(0, 2 * np.pi, 361)

        def pts(groups, color):
            out = []
            for c, m in groups:
                item = {"value": [r6(c.real), r6(c.imag)]}
                if m > 1:
                    item["label"] = {"show": True, "formatter": str(m), "position": "right", "color": color}
                out.append(item)
            return out

        cross = "path://M2,0 L5,3 L8,0 L10,2 L7,5 L10,8 L8,10 L5,7 L2,10 L0,8 L3,5 L0,2 Z"
        opt = {
            "animation": False,
            "grid": {"left": 60, "right": 30, "top": 40, "bottom": 50},
            "legend": {"top": 5},
            "tooltip": {"trigger": "item"},
            "toolbox": {"right": 10, "feature": {"dataZoom": {}, "restore": {}, "saveAsImage": {}}},
            "xAxis": {"type": "value", "min": -r6(rmax), "max": r6(rmax), "name": "Real",
                      "nameLocation": "middle", "nameGap": 28},
            "yAxis": {"type": "value", "min": -r6(rmax), "max": r6(rmax), "name": "Imaginary",
                      "nameLocation": "middle", "nameGap": 40},
            "series": [
                {"name": "Unit circle", "type": "line", "data": xy(np.cos(th), np.sin(th)), "showSymbol": False,
                 "silent": True, "lineStyle": {"color": P()["circle"], "type": "dashed", "width": 1},
                 "itemStyle": {"color": P()["circle"]}},
                {"name": "Zeros", "type": "scatter", "data": pts(cluster(z), P()["zero"]), "symbol": "circle",
                 "symbolSize": 11, "itemStyle": {"color": "rgba(0,0,0,0)", "borderColor": P()["zero"],
                                                 "borderWidth": 2}},
                {"name": "Poles", "type": "scatter", "data": pts(cluster(p), P()["pole"]), "symbol": cross,
                 "symbolSize": 11, "itemStyle": {"color": P()["pole"]}},
            ],
        }
        show(pz_chart, opt)

        pz_table.clear()
        with pz_table:
            ui.label("Poles").classes("font-semibold")
            for c, m in cluster(p):
                if c.imag >= -1e-12:
                    ui.label(f"r = {abs(c):.6f}   f = {abs(np.angle(c)) * fs / (2 * np.pi):9.2f} Hz"
                             + (f"   ×{m}" if m > 1 else "")).classes("font-mono")
            ui.label("Zeros").classes("font-semibold mt-2")
            for c, m in cluster(z):
                if c.imag >= -1e-12:
                    ui.label(f"r = {abs(c):.6f}   f = {abs(np.angle(c)) * fs / (2 * np.pi):9.2f} Hz"
                             + (f"   ×{m}" if m > 1 else "")).classes("font-mono")

    def render_time(d, impl_sos, fixed_imp, n) -> None:
        step = time_kind.value == "step"
        k = np.arange(n)
        y = fdcore.step_response(d.sos, n) if step else fdcore.impulse_response(d.sos, n)
        opt = base_chart("Sample n", "Amplitude")
        opt["series"].append(line("Design", xy(k, y), P()["design"]))
        if fixed_imp is not None:
            yi = fixed_imp
            opt["series"].append(line("Fixed point", xy(k, yi), P()["impl"], dashed=True, width=1.5))
        else:
            yi = fdcore.step_response(impl_sos, n) if step else fdcore.impulse_response(impl_sos, n)
            opt["series"].append(line("Implemented", xy(k, yi), P()["impl"], dashed=True, width=1.5))
        show(time_chart, opt)

    def render_sections(d, sos, fixed, fs) -> None:
        peaks = fdcore.partial_peak_gains(sos, fs)
        x = np.zeros(8192)
        x[0] = 1
        cols = [{"name": "k", "label": "#", "field": "k", "align": "right"},
                {"name": "ord", "label": "Order", "field": "ord", "align": "right"}]
        for nme in ("b0", "b1", "b2", "a1", "a2"):
            cols.append({"name": nme, "label": nme, "field": nme, "align": "right", "classes": "font-mono"})
        cols += [{"name": "pr", "label": "Pole r", "field": "pr", "align": "right"},
                 {"name": "pf", "label": "Pole f (Hz)", "field": "pf", "align": "right"},
                 {"name": "q", "label": "Pole Q", "field": "q", "align": "right"},
                 {"name": "pk", "label": "Peak gain at output", "field": "pk", "align": "right"},
                 {"name": "l1", "label": "L1 norm at output", "field": "l1", "align": "right"}]
        rows = []
        y = x
        for i, row in enumerate(sos):
            y = signal.sosfilt(row[None, :], y)
            p = np.roots(row[3:]) if row[5] != 0 else np.roots(row[3:5])
            pc = p[np.argmax(np.abs(p.imag))] if len(p) else 0
            r = abs(pc)
            theta = abs(np.angle(pc))
            # Q of the equivalent analog pole pair s = fs * ln(z)
            q = (math.hypot(math.log(r), theta) / (2 * abs(math.log(r)))
                 if row[5] != 0 and 0 < r < 1 and theta > 0 else None)
            if fixed is not None:
                vals = [str(int(v)) for v in fixed.coeffs[i]]
            else:
                vals = [f"{v:.9g}" for v in (row[0], row[1], row[2], row[4], row[5])]
            rows.append({"k": i + 1, "ord": d.section_order[i], **dict(zip(("b0", "b1", "b2", "a1", "a2"), vals)),
                         "pr": f"{r:.5f}", "pf": f"{theta * fs / (2 * np.pi):.1f}",
                         "q": f"{q:.2f}" if q else "–",
                         "pk": f"{20 * math.log10(peaks[i]):+.2f} dB",
                         "l1": f"{np.sum(np.abs(y)):.3f}"})
        sec_table.columns = cols
        sec_table.rows = rows
        sec_table.update()
        l1max = max(float(r["l1"]) for r in rows)
        sec_info.text = (
            ("Coefficients as integers in " + f"Q{fixed.frac} (int{fixed.word}_t). " if fixed is not None
             else "Coefficients in floating point. ")
            + ("Section scaling is on: the peak gain at every section output is 0 dB. " if scale_chk.value
               else "Section scaling is off. ")
            + f"Worst case input-to-node gain (L1 norm) is {l1max:.2f}, so "
            + f"{max(0, math.ceil(math.log2(l1max))) if l1max > 0 else 0} bit(s) of headroom avoid any overflow.")

    # ------------------------------------------------------------------ wiring
    def on_type_change(_=None):
        load_edges()
        schedule()

    def toggle_theme():
        dark.value = not dark.value
        theme_btn.props(f'icon={"light_mode" if dark.value else "dark_mode"}')
        update()

    theme_btn.on_click(toggle_theme)
    theme_btn.props(f'icon={"light_mode" if dark.value else "dark_mode"}')
    type_sel.on_value_change(on_type_change)
    frac_in.on_value_change(lambda: None if auto_q.value else schedule())
    for w in (char_sel, fs_in, ap_in, as_in, arith_sel, scale_chk, auto_q, name_in,
              mag_scale, mag_xlog, mag_mask, mag_impl, phase_wrap, gd_unit, time_kind, *edge_in):
        w.on_value_change(schedule)
    frac_in.bind_enabled_from(auto_q, "value", backward=lambda v: not v)

    load_edges()
    update()


def main() -> None:
    ap = argparse.ArgumentParser(description="FilterDesign GUI")
    ap.add_argument("--host", default="127.0.0.1", help="use 0.0.0.0 to allow access from other machines")
    ap.add_argument("--port", type=int, default=8090)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--light", action="store_true", help="start in light mode (default: dark)")
    args = ap.parse_args()
    global START_DARK
    START_DARK = not args.light
    ui.run(title="FilterDesign", host=args.host, port=args.port, show=not args.no_browser, reload=False, favicon="〰️")


if __name__ in {"__main__", "__mp_main__"}:
    main()
