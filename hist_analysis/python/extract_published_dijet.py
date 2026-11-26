"""Reconstruct the five legacy pPb 8.16 TeV ROOT-generated PDFs.

Run from the repository root with the project Python, CERN PyROOT and mplhep:
    PYTHONPATH="$(root-config --libdir)" py-env/bin/python -m hist_analysis.python.extract_published_dijet
"""

from array import array
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
import zlib

import numpy as np
import ROOT

from hist_analysis.python.plot_labels import PPB_MPLHEP_HEADER_RIGHT

BASE = Path(__file__).resolve().parents[1]
SOURCES = BASE / "root_files/old_distributions"
OUTPUT = SOURCES / "pPb8160_dijet_published_distributions.root"
PLOTS = BASE / "output/published_dijet_validation"
TOTAL = "pPb8160_syst_pt_all_lab.pdf"
FILES = {
    f"pPb8160_etaDijet_{tag}_h{half}_{frame}.pdf": observable
    for tag, frame, observable in [
        ("all", "lab", "full_eta"), ("fb", "cms", "forward_backward")
    ]
    for half in (1, 2)
}
NUM = r"-?\d+(?:\.\d+)?"
PAIR = rf"({NUM}) ({NUM})"
RECT = " ".join([f"({NUM})"] * 4)
FRAME = rf"0 0 0 RG 0 0 0 rg {RECT} re S "
SEGMENT = rf"{PAIR} m {PAIR} l S "
BAND = rf"0\.69 0\.81 0\.78 rg /ca0\.35 gs /CA0\.35 gs {RECT} re f\*"
# ROOT uses two cubic halves, with their shared extrema giving the center.
CIRCLE = rf"{PAIR} m (?:{NUM} ){{4}}{PAIR} c (?:{NUM} ){{6}}c h f"
MARKER = rf"1\.5 w ((?:{SEGMENT})*)0\.75 w {CIRCLE}"
TOL = 0.003  # PDF points, covers rounding of the operands and summed rectangle edges
AXES = {
    "full_eta": ([-2, 0, 2], [0.05, 0.1], "#eta^{dijet}_{lab}", "dN/d#eta"),
    "forward_backward": ([0, .5, 1, 1.5, 2], [1, 1.5, 2],
                         "|#eta^{dijet}_{CM}|", "Forward / Backward"),
    "total": ([-3, -2, -1, 0, 1, 2, 3], [0, 5, 10, 15, 20],
              "#eta^{dijet}_{lab}", "Relative total systematic uncertainty [%]"),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def page_stream(path):
    """Read the untransformed pad operands from the unique ROOT page stream."""
    streams = []
    for candidate in re.findall(rb"\bstream\r?\n(.*?)\r?\nendstream", path.read_bytes(), re.S):
        try:
            decoded = zlib.decompress(candidate).decode("latin1")
        except zlib.error:
            continue
        if re.search(FRAME, decoded):
            streams.append(decoded)
    require(len(streams) == 1, f"{path}: expected one ROOT graphics stream")
    require(streams[0].startswith("1 0 0 1 20 28 cm "), f"{path}: unexpected page transform")
    return streams[0]


def calibrate(labels, ticks):
    require(len(labels) == len(ticks), f"Unexpected tick sequence {ticks} for {labels}")
    slope, offset = np.polyfit(ticks, labels, 1)
    residual = max(abs(slope * np.array(ticks) + offset - labels))
    require(residual < abs(slope) * TOL, "Tick coordinates are inconsistent with a linear axis")
    return float(slope), float(offset)


@dataclass
class Panel:
    filename: str
    index: int
    pt: tuple
    kind: str
    frame: tuple
    xscale: tuple
    yscale: tuple
    content: str

    def x(self, value):
        return self.xscale[0] * value + self.xscale[1]

    def y(self, value):
        return self.yscale[0] * value + self.yscale[1]

    def edges(self, left, right):
        result = [round(self.x(edge), 1) for edge in (left, right)]
        require(max(abs(result[i] - self.x(edge)) for i, edge in enumerate((left, right)))
                < abs(self.xscale[0]) * TOL, f"Unexpected eta edges in {self.filename}, {self.pt}")
        return result


def panels(filename, kind):
    pieces = re.split(FRAME, page_stream(SOURCES / filename))
    expected = 20 if kind == "total" else 8
    require(len(pieces) == 1 + 5 * expected, f"{filename}: unexpected panel count")
    xlabels, ylabels, *_ = AXES[kind]
    result = []
    for index in range(expected):
        start = 1 + 5 * index
        left, bottom, width, height = map(float, pieces[start:start + 4])
        content = pieces[start + 4]
        low, = re.findall(r"\((\d+) < p\) Tj", content)
        high, = re.findall(r"\\\(GeV\\\) < (\d+)\) Tj", content)
        # Keep the literal precision of the PDF operands for the axis anchor.
        lx, by = map(re.escape, pieces[start:start + 2])
        horizontal_axis = re.search(rf"{lx} {by} m {NUM} {by} l S ", content)
        vertical_axis = re.search(rf"{lx} {by} m {lx} {NUM} l S ", content)
        require(horizontal_axis and vertical_axis, f"{filename}: axes missing")
        xticks = [(x, y) for x, y, xx, yy in segments(
            content[horizontal_axis.end():vertical_axis.start()])
            if x == xx and yy == bottom and y > bottom]
        # Minor ticks in the Total figure have half the length of major ticks.
        max_length = max(y - bottom for _, y in xticks)
        xs = sorted({x for x, y in xticks if y - bottom > 0.8 * max_length})
        yticks = segments(content[vertical_axis.end():].split(" q ", 1)[0])
        ys = sorted({y for x, y, xx, yy in yticks if x > left and xx == left and y == yy})
        result.append(Panel(filename, index + 1, (int(low), int(high)), kind,
                            (left, bottom, width, height), calibrate(xlabels, xs),
                            calibrate(ylabels, ys), content))
    return result


def segments(content):
    return [tuple(map(float, match)) for match in re.findall(SEGMENT, content + " ")]


def total_rows(panel):
    curve, = re.findall(rf"^1\.5 w ({PAIR} m (?:{PAIR} l )+)S ", panel.content)
    points = [tuple(map(float, v)) for v in re.findall(PAIR, curve[0])]
    rows = []
    for (x, y), (xx, yy) in zip(points, points[1:]):
        if xx > x:
            require(y == yy, f"Nonhorizontal Total step in {panel.pt}")
            percent = panel.y(y)
            clipped = abs(y - panel.frame[1] - panel.frame[3]) < TOL
            rows.append({"edges": panel.edges(x, xx), "percent": None if clipped else percent,
                         "drawn_percent": percent, "clipped_at_percent": 20 if clipped else None})
    require(len(rows) == 26 and all(r["drawn_percent"] > 0 for r in rows),
            f"Unexpected Total steps in {panel.pt}")
    require(all(a["edges"][1] == b["edges"][0] for a, b in zip(rows, rows[1:])),
            f"Duplicate or noncontiguous Total bins in {panel.pt}")
    return rows


def measurement_rows(panel):
    _, bottom, _, height = panel.frame
    rectangles = [tuple(map(float, b)) for b in re.findall(BAND, panel.content)]
    data_section = panel.content.split(" Tj ET Q 1.5 w ", 1)[1]
    data_section = "1.5 w " + data_section.split(" q ", 1)[0]
    markers = list(re.finditer(MARKER, data_section))
    require(len(markers) == len(rectangles), f"{panel.filename} {panel.pt}: band/marker count")
    rows = []
    for marker, (bx, by, bw, bh) in zip(markers, rectangles):
        px, py, qx, qy = map(float, marker.groups()[-4:])
        cx = (px + qx) / 2
        require(abs(py - qy) < TOL and abs(cx - bx - bw / 2) < TOL,
                f"{panel.filename} {panel.pt}: band/marker alignment")
        edges = panel.edges(bx, bx + bw)
        y = panel.y(py)
        x = sum(edges) / 2
        require(abs(panel.x(cx) - x) < abs(panel.xscale[0]) * TOL,
                f"{panel.filename} {panel.pt}: eta midpoint mismatch")
        endpoints = [v for x1, y1, x2, y2 in segments(marker[1])
                     if abs(x1 - cx) < TOL and abs(x2 - cx) < TOL for v in (y1, y2)]
        require(len(endpoints) in (0, 2, 4), f"Unexpected statistical path in {panel.pt}")
        for a, b in zip(endpoints[::2], endpoints[1::2]):
            require(min(abs(a - py), abs(b - py)) <= 3.6306 + TOL
                    and (a - py) * (b - py) > 0, f"Unexpected marker gap in {panel.pt}")
        stat = []
        bounds = []
        stat_drawn = []
        for sign, boundary in [(-1, bottom), (1, bottom + height)]:
            distances = [sign * (endpoint - py) for endpoint in endpoints
                         if sign * (endpoint - py) > 0]
            distance = max(distances, default=None)
            clipped = distance is not None and abs(py + sign * distance - boundary) < TOL
            stat_drawn.append(None if distance is None else distance * panel.yscale[0])
            stat.append(None if distance is None or clipped else distance * panel.yscale[0])
            # A hidden side extends at most as far as the ROOT marker gap,
            # provided the frame leaves enough room to see the endpoint.
            gap = 3.6306
            bounds.append((gap + TOL) * panel.yscale[0] if distance is None
                          and sign * (boundary - py) > gap + TOL else None)
        syst = [py - by, by + bh - py]
        # An edge on the frame can hide one side of a symmetric rectangle.
        clipped_sides = [abs(by - bottom) < TOL, abs(by + bh - bottom - height) < TOL]
        if any(clipped_sides):
            require(clipped_sides.count(True) == 1, f"Both band edges clipped in {panel.pt}")
            visible = syst[clipped_sides.index(False)]
            require(visible + TOL >= syst[clipped_sides.index(True)],
                    f"Band clipping incompatible with symmetry in {panel.pt}")
            syst = [visible, visible]
        require(min(syst) >= 0 and abs(syst[0] - syst[1]) < TOL,
                f"Band is inconsistent with symmetry in {panel.pt}, eta={x}")
        band_error = sum(syst) / 2 * panel.yscale[0]
        require(y > 0, f"Expected positive central value in {panel.pt}")
        rows.append({"edges": edges, "x": x, "y": y, "stat": stat,
                     "stat_hidden_upper": bounds, "stat_drawn": stat_drawn,
                     "band_error": band_error,
                     "band_clipped": any(clipped_sides),
                     "syst": band_error, "relative": band_error / y,
                     "syst_source": "green measurement band"})
    require(all(a["edges"][1] == b["edges"][0] for a, b in zip(rows, rows[1:])),
            f"Duplicate or noncontiguous eta bins in {panel.pt}")
    return rows


def collect():
    totals = {panel.pt: (panel, total_rows(panel)) for panel in panels(TOTAL, "total")}
    measurements = [(panel, measurement_rows(panel))
                    for name, kind in FILES.items() for panel in panels(name, kind)]
    for kind in dict.fromkeys(FILES.values()):
        intervals = [p.pt for p, _ in measurements if p.kind == kind]
        require(len(intervals) == 16 and all(a[1] == b[0] for a, b in zip(intervals, intervals[1:]))
                and intervals[0][0] == 50 and intervals[-1][1] == 500,
                f"Unexpected pT coverage for {kind}")
    for panel, rows in measurements:
        if panel.kind != "full_eta" or panel.pt not in totals:
            continue
        _, total = totals[panel.pt]
        by_edges = {tuple(row["edges"]): row["percent"] for row in total}
        for row in rows:
            percent = by_edges[tuple(row["edges"])]
            require(percent is not None, f"Clipped matched Total in {panel.pt}")
            row["black_total_relative"] = percent / 100
            row["black_total_absolute"] = abs(row["y"]) * percent / 100
    return measurements, list(totals.values())


def directory_name(panel):
    base = "standalone_full_eta_total" if panel.kind == "total" else panel.kind
    return f"{base}/ptAve_{panel.pt[0]}_{panel.pt[1]}"


def graph_specs(panel, rows):
    # Tuples are x, y, horizontal low/high and vertical low/high errors.
    def point(row, value, low=0., high=0.):
        x = sum(row["edges"]) / 2
        return (x, value, x - row["edges"][0], row["edges"][1] - x, low, high)

    if panel.kind == "total":
        return {
            "data_syst_relative": [point(r, r["percent"] / 100) for r in rows
                                   if r["percent"] is not None],
            "displayed_total_percent": [point(r, r["drawn_percent"]) for r in rows],
        }
    result = {
        "data": [point(r, r["y"]) for r in rows],
        "data_stat": [point(r, r["y"], *r["stat"]) for r in rows if None not in r["stat"]],
        "data_syst": [point(r, r["y"], r["syst"], r["syst"]) for r in rows],
        "data_syst_relative": [point(r, r["relative"]) for r in rows],
    }
    if "black_total_relative" in rows[0]:
        result["black_total_relative"] = [point(r, r["black_total_relative"]) for r in rows]
        result["black_total_absolute"] = [point(r, r["y"], r["black_total_absolute"], r["black_total_absolute"])
                                          for r in rows]
    return result


def json_object(value):
    return ROOT.TObjString(json.dumps(value, allow_nan=False, sort_keys=True))


def save_root(path, records):
    with ROOT.TFile.Open(str(path), "RECREATE") as output:
        output.mkdir("metadata").WriteObject(json_object({
            "source_sha256": {name: sha256((SOURCES / name).read_bytes()).hexdigest()
                              for name in (*FILES, TOTAL)},
            "method": "ROOT PDF vector operands and linear major-tick calibration",
            "source": "CMS-PAS-HIN-24-014, measurement Figures 1–4 and local Total figure",
            "frame": {"full_eta": "beam-aligned laboratory eta",
                      "forward_backward": "positive (proton-going) CM eta / reflected negative CM eta"},
            "data": "all displayed central values",
            "data_stat": "only points with both statistical endpoints resolved; see per-panel rows",
            "stat_sides": ["lower", "upper"],
            "data_syst": "same-panel symmetric green systematic band",
            "data_syst_relative": "dimensionless fraction",
            "black_total": "alternative uncertainty from exact pT and eta matches; differs from green bands",
            "standalone_total": "displayed_total_percent includes steps clipped at 20%; "
                                "data_syst_relative excludes clipped values",
            "digitization": "coordinate tolerance of 0.003 PDF points, separate from physics errors",
        }), "provenance")
        for panel, rows in records:
            directory = output.mkdir(directory_name(panel))
            directory.WriteObject(json_object({
                "file": panel.filename, "page": 1, "panel": panel.index, "pt_gev": panel.pt,
                "pdf_frame": panel.frame, "x_calibration": panel.xscale, "y_calibration": panel.yscale,
                "rows": rows,
            }), "provenance")
            for name, points in graph_specs(panel, rows).items():
                graph = ROOT.TGraphAsymmErrors(len(points))
                graph.SetName(name)
                _, _, xlabel, ylabel = AXES[panel.kind]
                if name.endswith("relative"):
                    ylabel = "Relative total systematic uncertainty (fraction)"
                graph.SetTitle(f"{name}, {panel.pt[0]}–{panel.pt[1]} GeV;{xlabel};{ylabel}")
                for i, (x, y, xl, xh, yl, yh) in enumerate(points):
                    graph.SetPoint(i, x, y)
                    graph.SetPointError(i, xl, xh, yl, yh)
                directory.WriteObject(graph, name)


def read_points(graph):
    require(graph and graph.InheritsFrom("TGraphAsymmErrors"), "Missing ROOT graph")
    getters = (graph.GetPointX, graph.GetPointY, graph.GetErrorXlow, graph.GetErrorXhigh,
               graph.GetErrorYlow, graph.GetErrorYhigh)
    return [tuple(get(i) for get in getters) for i in range(graph.GetN())]


def validate_root(archive, records):
    provenance = json.loads(str(archive.Get("metadata/provenance").GetString()))
    for name, digest in provenance["source_sha256"].items():
        require(sha256((SOURCES / name).read_bytes()).hexdigest() == digest, f"Source changed: {name}")
    for panel, rows in records:
        directory = archive.GetDirectory(directory_name(panel))
        stored = json.loads(str(directory.Get("provenance").GetString()))
        require(stored["rows"] == rows and stored["pt_gev"] == list(panel.pt),
                f"ROOT metadata round trip failed for {directory_name(panel)}")
        for name, expected in graph_specs(panel, rows).items():
            actual = read_points(directory.Get(name))
            require(actual == expected, f"ROOT numerical round trip failed for {directory_name(panel)}/{name}")
            require(all(min(p[2:]) >= 0 for p in actual), "Negative error extent")
        if panel.kind != "total":
            syst = read_points(directory.Get("data_syst"))
            relative = read_points(directory.Get("data_syst_relative"))
            require(all(np.isclose(abs(a[1]) * b[1], a[-1], rtol=1e-14, atol=0)
                        for a, b in zip(syst, relative)), "Inconsistent absolute/relative errors")


def draw_measurement(ax, directory):
    for x, y, xl, xh, yl, yh in read_points(directory.Get("data_syst")):
        ax.fill_between([x - xl, x + xh], y - yl, y + yh, color="#b0cfc7",
                        alpha=.45, linewidth=0)
    data = np.array(read_points(directory.Get("data")))
    ax.plot(data[:, 0], data[:, 1], "o", color="black", markersize=3)
    for x, y, xl, xh, *_ in data:
        ax.hlines(y, x - xl, x + xh, color="black", linewidth=.65)
    for x, y, _, _, yl, yh in read_points(directory.Get("data_stat")):
        ax.vlines(x, y - yl, y + yh, color="black", linewidth=.8)
    for row in json.loads(str(directory.Get("provenance").GetString()))["rows"]:
        if None not in row["stat"]:
            continue
        # Include visible portions of frame-clipped bars for source comparison.
        lower, upper = row["stat_drawn"]
        if lower is not None or upper is not None:
            ax.vlines(row["x"], row["y"] - (lower or 0), row["y"] + (upper or 0),
                      color="black", linewidth=.8)


def validation_plots(archive, records, source_renders):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep
    from PIL import Image

    PLOTS.mkdir(parents=True, exist_ok=True)
    for filename in (*FILES, TOTAL):
        group = [panel for panel, _ in records if panel.filename == filename]
        is_total = filename == TOTAL
        nrows, ncols = (5, 4) if is_total else (4, 2)
        with plt.style.context(hep.style.CMS):
            fig, axes = plt.subplots(nrows, ncols, figsize=(12, 11), layout="constrained")
            fig.supxlabel(r"$\eta^{\mathrm{dijet}}_{\mathrm{lab}}$" if group[0].kind != "forward_backward"
                          else r"$|\eta^{\mathrm{dijet}}_{\mathrm{CM}}|$", fontsize=18)
            fig.supylabel({"total": "Relative uncertainty [%]", "full_eta": r"$\mathrm{d}N/\mathrm{d}\eta$",
                           "forward_backward": r"Forward ($\eta_{\rm CM}>0$) / Backward ($\eta_{\rm CM}<0$)"
                           }[group[0].kind], fontsize=16)
            if not is_total:
                hep.cms.label("Preliminary", data=True, rlabel=PPB_MPLHEP_HEADER_RIGHT,
                              ax=axes[0, 0], fontsize=14)
            for ax, panel in zip(axes.flat, group):
                directory = archive.GetDirectory(directory_name(panel))
                if is_total:
                    points = read_points(directory.Get("displayed_total_percent"))
                    edges = [p[0] - p[2] for p in points] + [points[-1][0] + points[-1][3]]
                    ax.stairs([p[1] for p in points], edges, color="black", linewidth=1.25)
                else:
                    draw_measurement(ax, directory)
                left, bottom, width, height = panel.frame
                ax.set(xlim=(panel.x(left), panel.x(left + width)),
                       ylim=(panel.y(bottom), panel.y(bottom + height)))
                ax.set_xticks(AXES[panel.kind][0])
                ax.set_yticks(AXES[panel.kind][1])
                ax.tick_params(labelsize=11)
                ax.text(.5, .95, rf"${panel.pt[0]} < p_{{\mathrm{{T}}}}^{{\mathrm{{ave}}}}"
                        rf" < {panel.pt[1]}$ GeV", transform=ax.transAxes, ha="center", va="top", fontsize=12)
            target = PLOTS / f"{Path(filename).stem}_comparison.png"
            reconstructed = source_renders / "reconstructed.png"
            fig.savefig(reconstructed, dpi=120)
            plt.close(fig)
        with Image.open(source_renders / f"{filename}.png") as source, Image.open(reconstructed) as redrawn:
            original = source.resize((round(source.width * redrawn.height / source.height), redrawn.height))
            comparison = Image.new("RGB", (original.width + redrawn.width, redrawn.height), "white")
            comparison.paste(original, (0, 0))
            comparison.paste(redrawn, (original.width, 0))
            comparison.save(target)


def raster_check(records, source_renders):
    """Independently test vector-derived marker centers against Poppler pixels."""
    from PIL import Image

    for panel, rows in records:
        if panel.kind == "total":
            continue
        with Image.open(source_renders / f"{panel.filename}.png") as source:
            pixels = np.asarray(source.convert("RGB"))
        # Render at 144 dpi: 2 px/PDF pt. PDF stream translates by (20, 28).
        for row in rows:
            x = round(2 * ((row["x"] - panel.xscale[1]) / panel.xscale[0] + 20))
            y = round(pixels.shape[0] - 2 * ((row["y"] - panel.yscale[1]) / panel.yscale[0] + 28))
            require(np.max(pixels[y - 1:y + 2, x - 1:x + 2]) < 40,
                    f"Poppler marker cross-check failed for {panel.filename} {panel.pt}, eta={row['x']}")


def main():
    ROOT.gROOT.SetBatch(True)
    measurements, totals = collect()
    records = measurements + totals
    with TemporaryDirectory(dir=SOURCES) as temporary:
        temp = Path(temporary)
        candidate = temp / OUTPUT.name
        save_root(candidate, records)
        with ROOT.TFile.Open(str(candidate), "READ") as archive:
            validate_root(archive, records)
            for filename in (*FILES, TOTAL):
                subprocess.run(["pdftoppm", "-r", "144", "-singlefile", "-png",
                                str(SOURCES / filename), str(temp / filename)], check=True)
            raster_check(records, temp)
            validation_plots(archive, records, temp)
        candidate.replace(OUTPUT)
    print(f"Verified {len(measurements)} measurement panels and {len(totals)} Total panels in {OUTPUT}")


if __name__ == "__main__":
    main()
