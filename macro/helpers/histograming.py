import numpy as np
import ROOT
from ROOT import TH2D, TH1D, TCanvas, gStyle, THnSparseD
import array
import uproot
import sys, os
import getpass
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
import ctypes

from mpl_toolkits.axes_grid1 import make_axes_locatable
from pathlib import Path

# import mplhep as hep
# hep.style.use("CMS")

# Set global matplotlib style
rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.labelweight": "normal",
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "mathtext.default": "regular",  # use regular math font
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
})

# Helper function to set default style
def setDefaultStyle():
    """Configure base ROOT style."""
    gStyle = ROOT.gStyle
    if not gStyle:
        return
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetPalette(ROOT.kBird)

# Helper function to detect trigger from filename
def detect_trigger(filename):
    """Infer trigger type from filename."""
    if "Jet100" in filename:
        return "Jet100"
    elif "Jet80" in filename:
        return "Jet80"
    elif "Jet60" in filename:
        return "Jet60"
    elif "MB" in filename:
        return "MB"
    else:
        return "Unknown"

# Helper function to detect direction from filename
def detect_direction(filename):
    """Infer direction (Pbgoing, pgoing, combined) from filename."""
    if "Pbgoing" in filename:
        return "Pbgoing"
    elif "pgoing" in filename:
        return "pgoing"
    else:
        return "combined"

# Helper function to draw CMS header
def plotCMSHeader(collSystem=0, energy=5.02):
    """Draw CMS header. collSystem: 0=pp, 1=pPb, 2=PbPb; energy in TeV"""
    collSystemStr = "pp" if collSystem == 0 else ("pPb" if collSystem == 1 else "PbPb")
    t = ROOT.TLatex()
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.DrawLatexNDC(0.15, 0.93, "#bf{CMS} #it{Preliminary}")
    t.SetTextSize(0.04)
    t.DrawLatexNDC(0.6, 0.93, f"{collSystemStr} #sqrt{{s_{{NN}}}} = {energy:3.2f} TeV")
    t.SetTextSize(0.05)

# Helper function to set 1D histogram style
def set1DStyle(h, style_type=0, doRenorm=False):
    """Apply 1D styling to a ROOT.TH1"""
    markerStyle = 20
    markerSize = 1.2
    lineWidth = 2
    # default color = red
    color = ROOT.kRed

    if style_type == 0:
        color = ROOT.kRed
        markerStyle = 20
    elif style_type == 1:
        color = ROOT.kBlue
        markerStyle = 24
    elif style_type == 2:
        color = ROOT.kBlack
        markerStyle = 20
    elif style_type == 3:
        color = ROOT.kRed
        markerStyle = 24
    elif style_type == 4:
        color = ROOT.kBlue
        markerStyle = 20
    else:
        color = ROOT.kMagenta
        markerStyle = 30

    h.SetLineWidth(lineWidth)
    h.SetLineColor(color)

    h.SetMarkerStyle(markerStyle)
    h.SetMarkerColor(color)
    h.SetMarkerSize(markerSize)

    yaxis = h.GetYaxis()
    xaxis = h.GetXaxis()

    yaxis.SetTitleSize(0.05)
    yaxis.SetLabelSize(0.05)
    xaxis.SetTitleSize(0.05)
    xaxis.SetLabelSize(0.05)
    xaxis.SetNdivisions(208)
    yaxis.SetNdivisions(208)
    yaxis.SetTitleOffset(1.1)

    if doRenorm:
        integral = h.Integral()
        if integral != 0:
            h.Scale(1.0 / integral)

# Helper function to set 2D histogram style
def set2DStyle(h):
    """Apply 2D styling to a ROOT.TH2"""
    yaxis = h.GetYaxis()
    xaxis = h.GetXaxis()

    yaxis.SetTitleSize(0.06)
    yaxis.SetLabelSize(0.06)
    xaxis.SetTitleSize(0.06)
    xaxis.SetLabelSize(0.06)
    xaxis.SetNdivisions(205)
    yaxis.SetNdivisions(205)
    yaxis.SetTitleOffset(1.0)

# Helper function to set standard pad margins
def setPadStyle():
    """Set standard pad margins on the current ROOT pad"""
    pad = ROOT.gPad
    if not pad:
        return
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.15)
    pad.SetRightMargin(0.1)
    pad.SetLeftMargin(0.15)
    pad.SetGrid()

# Helper function to rescale 1D histogram
def rescaleHisto1D(h1):
    """Rescale a 1D ROOT.TH1 so the area (sum content*binWidth) is normalized and then enforce unit integral."""
    if not h1:
        sys.stderr.write("Null histogram pointer!\n")
        return

    nBins = h1.GetNbinsX()
    xaxis = h1.GetXaxis()

    total = 0.0
    for i in range(1, nBins + 1):
        width = xaxis.GetBinWidth(i)
        content = h1.GetBinContent(i)
        total += content * width

    if total > 0.0:
        for i in range(1, nBins + 1):
            content = h1.GetBinContent(i)
            error = h1.GetBinError(i)
            h1.SetBinContent(i, content / total)
            h1.SetBinError(i, error / total)
    else:
        sys.stderr.write("Warning: total area-normalized sum is zero!\n")

    integral = h1.Integral()
    if integral != 0.0:
        h1.Scale(1.0 / integral)

# Helper function to add CMS header to matplotlib axes
def add_cms_header(ax, label_type="Preliminary", energy_label=r"pPb $\sqrt{s_{\mathrm{NN}}}$ = 8.16 TeV"):
    """
    Add CMS-style header to the top of a plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to annotate.
    label_type : str
        'Preliminary', 'Simulation', or custom text. Italicized next to CMS.
    energy_label : str
        Right-side text, e.g. collision system and energy.
    """

    # --- Build left header text (CMS bold, suffix italic) ---
    # Use mathtext formatting to mix bold and italic in the same line
    cms_text = r"$\mathbf{CMS}$"
    if label_type:
        cms_text += rf" $\it{{{label_type}}}$"

    # --- Compute the position of the axes in figure coordinates ---
    pos = ax.get_position()
    y_offset = 0.015  # vertical offset above plot box

    # --- Left text: 'CMS Preliminary' ---
    ax.figure.text(
        pos.x0, pos.y1 + y_offset,
        cms_text,
        ha="left", va="bottom",
        fontsize=16,
        fontfamily="sans-serif",
        usetex=False
    )

    # --- Right text: 'pPb √sₙₙ = 8.16 TeV' ---
    ax.figure.text(
        pos.x1, pos.y1 + y_offset,
        energy_label,
        ha="right", va="bottom",
        fontsize=16,
        fontfamily="sans-serif",
        usetex=False
    )

# Helper function to read histogram from ROOT file
def read_histogram(path, hist_name):
    """
    Reads a histogram (1D, 2D, 3D, or THnSparseD) from a ROOT file.
    Works with both uproot and ROOT.
    Returns a consistent dictionary for downstream analysis.

    Returns
    -------
    dict with keys:
        ndim: int
        values: np.ndarray
        errors: np.ndarray
        edges: list of np.ndarray  (each axis bin edges)
        labels: list of str        (axis titles)
    """

    # --- Try uproot first (fast + safe) ---
    try:
        with uproot.open(path) as f:
            if hist_name not in f:
                raise KeyError(f"Histogram '{hist_name}' not found in file {path}")
            h = f[hist_name]

            ndim = len(h.axes)
            if ndim == 0:
                raise TypeError(f"Unsupported histogram structure: {type(h)}")

            # --- Retrieve values and variances ---
            values = h.values(flow=False)
            variances = h.variances(flow=False)
            errors = np.sqrt(variances) if variances is not None else np.zeros_like(values)

            # --- Axis bin edges and labels ---
            edges = [ax.edges() for ax in h.axes]
            labels = [ax.title or f"axis_{i}" for i, ax in enumerate(h.axes)]

            return {
                "ndim": ndim,
                "values": values,
                "errors": errors,
                "edges": edges,
                "labels": labels,
                "source": "uproot",
            }

    except Exception as e:
        print(f"[uproot] Failed to read '{hist_name}' from {path}: {e}")
        # Fallback to PyROOT (for THnSparse, 3D, or special cases)

    # --- ROOT fallback ---
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise IOError(f"Cannot open ROOT file: {path}")

    h = f.Get(hist_name)
    if not h:
        f.Close()
        raise KeyError(f"Histogram '{hist_name}' not found in {path}")

    # --- TH1D / TH2D / TH3D ---
    if h.InheritsFrom("TH1"):
        nx = h.GetNbinsX()
        values = np.array([h.GetBinContent(i) for i in range(1, nx + 1)])
        errors = np.array([h.GetBinError(i) for i in range(1, nx + 1)])
        edges = [np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])]
        labels = [h.GetXaxis().GetTitle() or "X"]

        result = {"ndim": 1, "values": values, "errors": errors,
                  "edges": edges, "labels": labels, "source": "ROOT"}

    elif h.InheritsFrom("TH2"):
        nx, ny = h.GetNbinsX(), h.GetNbinsY()
        values = np.zeros((ny, nx))
        errors = np.zeros_like(values)
        for ix in range(1, nx + 1):
            for iy in range(1, ny + 1):
                values[iy - 1, ix - 1] = h.GetBinContent(ix, iy)
                errors[iy - 1, ix - 1] = h.GetBinError(ix, iy)
        xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
        yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
        result = {"ndim": 2, "values": values, "errors": errors,
                  "edges": [xedges, yedges],
                  "labels": [h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()],
                  "source": "ROOT"}

    elif h.InheritsFrom("TH3"):
        nx, ny, nz = h.GetNbinsX(), h.GetNbinsY(), h.GetNbinsZ()
        values = np.zeros((nz, ny, nx))
        errors = np.zeros_like(values)
        for ix in range(1, nx + 1):
            for iy in range(1, ny + 1):
                for iz in range(1, nz + 1):
                    values[iz - 1, iy - 1, ix - 1] = h.GetBinContent(ix, iy, iz)
                    errors[iz - 1, iy - 1, ix - 1] = h.GetBinError(ix, iy, iz)
        xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
        yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
        zedges = np.array([h.GetZaxis().GetBinLowEdge(i) for i in range(1, nz + 2)])
        result = {"ndim": 3, "values": values, "errors": errors,
                  "edges": [xedges, yedges, zedges],
                  "labels": [h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle(), h.GetZaxis().GetTitle()],
                  "source": "ROOT"}

    elif h.InheritsFrom("THnSparse"):
        ndim = h.GetNdimensions()
        shape = [h.GetAxis(i).GetNbins() for i in range(ndim)]
        values = np.zeros(shape)
        errors = np.zeros_like(values)

        # Loop over all filled bins
        it = h.GetSparseIterator()
        while True:
            bin_idx = it.Next()
            if bin_idx < 0:
                break
            content = h.GetBinContent(bin_idx)
            error = h.GetBinError(bin_idx)
            indices = [h.GetBinCoordinates(bin_idx)[i] - 1 for i in range(ndim)]
            try:
                values[tuple(indices)] = content
                errors[tuple(indices)] = error
            except Exception:
                pass

        edges = [np.array([h.GetAxis(i).GetBinLowEdge(j)
                           for j in range(1, h.GetAxis(i).GetNbins() + 2)]) for i in range(ndim)]
        labels = [h.GetAxis(i).GetTitle() or f"axis_{i}" for i in range(ndim)]

        result = {"ndim": ndim, "values": values, "errors": errors,
                  "edges": edges, "labels": labels, "source": "ROOT_THnSparse"}

    else:
        f.Close()
        raise TypeError(f"Unsupported object type: {type(h)}")

    f.Close()
    return result

# Helper function to plot 2D histogram with CMS-style formatting
def plot_2d_histogram(
    h,
    outname,
    xlabel=r"$x$",
    ylabel=r"$y$",
    show_plot=False,
    is_preliminary=True,
    xlim=None,
    ylim=None,
    energy_label=r"pPb $\sqrt{s_{\mathrm{NN}}}$ = 8.16 TeV",
    cmap="viridis",
    cbar_label=None
):
    """
    Plot any 2D histogram (ROOT or uproot) with CMS-style formatting.

    Parameters
    ----------
    h : uproot.behaviors.TH2.Histogram or ROOT.TH2
        Input 2D histogram.
    outname : str
        Output filename (e.g. "plots/h2D_example.pdf").
    xlabel : str
        X-axis label (LaTeX allowed).
    ylabel : str
        Y-axis label (LaTeX allowed).
    show_plot : bool
        If True, display plot inline; otherwise save and close.
    is_preliminary : bool
        Controls CMS header text ("Preliminary" vs "Simulation").
    xlim, ylim : tuple(float, float), optional
        Axis range limits.
    energy_label : str
        Text to appear in the top-right corner (default: pPb sqrt{s} = 8.16 TeV).
    cmap : str
        Matplotlib colormap name (default: "viridis").
    cbar_label : str
        Label for the colorbar.
    """

    # --- Convert ROOT/uproot histogram to numpy arrays ---
    if hasattr(h, "to_numpy"):
        vals, xedges, yedges = h.to_numpy()
    elif hasattr(h, "GetArray"):
        # Support for ROOT.TH2 histograms
        import numpy as np
        nx, ny = h.GetNbinsX(), h.GetNbinsY()
        vals = np.zeros((ny, nx))
        for ix in range(1, nx + 1):
            for iy in range(1, ny + 1):
                vals[iy - 1, ix - 1] = h.GetBinContent(ix, iy)
        xedges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
        yedges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
    else:
        raise TypeError("Unsupported histogram type — must be uproot or ROOT.TH2")

    # --- Create a square figure ---
    fig, ax = plt.subplots(figsize=(6, 6))

    # --- Draw the 2D histogram ---
    mesh = ax.pcolormesh(xedges, yedges, vals.T, cmap=cmap)

    # --- Axis labels ---
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    # --- Axis limits (if provided) ---
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    # --- Equal aspect ratio ---
    ax.set_box_aspect(1)

    # --- Colorbar aligned with histogram height ---
    pos = ax.get_position()
    cbar_width = 0.03
    cbar_padding = 0.01
    cax = fig.add_axes([
        pos.x1 + cbar_padding,
        pos.y0,
        cbar_width,
        pos.height
    ])
    cbar = fig.colorbar(mesh, cax=cax)
    if cbar_label:
        cbar.set_label(cbar_label, fontsize=14)

    # --- Axis styling ---
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=13)

    # --- Add CMS header ---
    label_type = "Preliminary" if is_preliminary else "Simulation"
    add_cms_header(ax, label_type=label_type, energy_label=energy_label)

    # --- Save or show ---
    fig.savefig(outname, bbox_inches="tight")
    if show_plot:
        plt.show()
    else:
        plt.close(fig)

# Helper function to plot 1D histogram with CMS-style formatting
def plot_1d_histogram(
    h,
    outname,
    xlabel=r"$x$",
    ylabel=r"Entries",
    show_plot=False,
    is_preliminary=True,
    xlim=None,
    ylim=None,
    energy_label=r"pPb $\sqrt{s_{\mathrm{NN}}}$ = 8.16 TeV",
    color="black",
    marker="o",
    markersize=4,
    linestyle="-",
    linewidth=1.5,
    grid=False,
):
    """
    Plot any 1D histogram (ROOT or uproot) with CMS-style formatting.
    Includes full bin error support for both ROOT.TH1 and uproot histograms.

    Parameters
    ----------
    h : uproot.behaviors.TH1.Histogram or ROOT.TH1
        Input 1D histogram.
    outname : str
        Output filename (e.g. "plots/h1D_example.pdf").
    xlabel, ylabel : str
        Axis labels (LaTeX allowed).
    show_plot : bool
        If True, show inline; otherwise save to file.
    is_preliminary : bool
        Controls CMS header ("Preliminary" vs "Simulation").
    xlim, ylim : tuple(float, float)
        Axis limits.
    energy_label : str
        Top-right text (collision system and energy).
    color, marker, linestyle : str
        Matplotlib style controls.
    grid : bool
        Enable grid lines.
    """

    # --- Convert histogram to numpy arrays ---
    # uproot histogram
    if hasattr(h, "to_numpy"):
        values, edges = h.to_numpy()
        centers = 0.5 * (edges[1:] + edges[:-1])
        # Get bin errors
        if hasattr(h, "variances") and h.variances() is not None:
            errors = np.sqrt(h.variances())
        else:
            errors = np.zeros_like(values)
    # ROOT histogram
    elif hasattr(h, "GetNbinsX"):
        nbins = h.GetNbinsX()
        values = np.array([h.GetBinContent(i) for i in range(1, nbins + 1)])
        centers = np.array([h.GetXaxis().GetBinCenter(i) for i in range(1, nbins + 1)])
        errors = np.array([h.GetBinError(i) for i in range(1, nbins + 1)])
    else:
        raise TypeError("Unsupported histogram type — must be uproot or ROOT.TH1")

    # --- Create figure ---
    fig, ax = plt.subplots(figsize=(6, 5))

    # --- Draw histogram with error bars ---
    ax.errorbar(
        centers,
        values,
        yerr=errors,
        fmt=marker,
        color=color,
        markersize=markersize,
        linestyle=linestyle,
        linewidth=linewidth,
        capsize=2,
        elinewidth=1.0,
    )

    # --- Axis labels ---
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)

    # --- Axis limits ---
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    # --- Axis formatting ---
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=13)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    if grid:
        ax.grid(True, linestyle="--", alpha=0.4)

    # --- Add CMS header ---
    label_type = "Preliminary" if is_preliminary else "Simulation"
    add_cms_header(ax, label_type=label_type, energy_label=energy_label)

    # --- Save or show ---
    fig.savefig(outname, bbox_inches="tight")
    if show_plot:
        plt.show()
    else:
        plt.close(fig)

# Helper function to project TH2 onto Y axis for given X interval
def projectDijetEtaFrom2D(h2, x_low=None, x_high=None):
    """
    Project a 2D ROOT.TH2 histogram onto the Y axis for a given X range.

    Parameters
    ----------
    h2 : ROOT.TH2
        Input 2D histogram.
    x_low : float or int, optional
        Lower edge or bin number along X-axis to start projection.
    x_high : float or int, optional
        Upper edge or bin number along X-axis to end projection.

    Returns
    -------
    ROOT.TH1D
        Projected histogram along Y axis.
    """
    # Check input
    if not h2 or not h2.InheritsFrom("TH2"):
        sys.stderr.write("Error: Input must be a valid TH2 histogram\n")
        return None

    xaxis = h2.GetXaxis()
    nxbins = xaxis.GetNbins()

    # Convert values to bin indices
    def to_bin(v, default_bin):
        if v is None:
            return default_bin
        if isinstance(v, int):
            return max(1, min(nxbins, v))
        try:
            return max(1, min(nxbins, xaxis.FindBin(float(v))))
        except Exception:
            return default_bin

    bin_low = to_bin(x_low, 1)
    bin_high = to_bin(x_high, nxbins)

    # Ensure proper ordering
    if bin_low > bin_high:
        bin_low, bin_high = bin_high, bin_low

    # Perform projection
    proj_name = f"{h2.GetName()}_ptAve_{bin_low}_{bin_high}"
    proj = h2.ProjectionY(proj_name, bin_low, bin_high)

    if not proj:
        sys.stderr.write("Error: Projection failed\n")
        return None

    proj.SetDirectory(0)  # detach from file
    proj.GetXaxis().SetTitle(h2.GetYaxis().GetTitle())
    proj.GetYaxis().SetTitle("Entries")

    return proj

# Helper function to create Forward/Backward ratio histogram
def make_fb_ratio(h_forward, h_backward):
    """
    Create a Forward/Backward ratio histogram.

    Parameters
    ----------
    h_forward : ROOT.TH1
        Histogram for the forward region. Its name should contain 'Forward'.
    h_backward : ROOT.TH1
        Histogram for the backward region.

    Returns
    -------
    ROOT.TH1
        Ratio histogram (Forward / Backward) with the name where 'Forward'
        is replaced by 'FBRatio'.
    """
    if not h_forward or not h_backward:
        raise ValueError("Both forward and backward histograms must be provided.")
    if not (h_forward.InheritsFrom("TH1") and h_backward.InheritsFrom("TH1")):
        raise TypeError("Inputs must be TH1 histograms.")

    # Clone forward histogram
    new_name = h_forward.GetName().replace("Forward", "FBRatio")
    h_ratio = h_forward.Clone(new_name)
    h_ratio.SetDirectory(0)  # prevent deletion when file closes

    # Compute ratio safely
    h_ratio.Divide(h_backward)

    # Set appearance
    h_ratio.SetTitle(h_forward.GetTitle().replace("Forward", "FBRatio"))
    set1DStyle(h_ratio, style_type=2)
    h_ratio.GetYaxis().SetTitle("Forward / Backward")

    return h_ratio

# Helper function to load histograms from ROOT file with tagging
def load_histograms(root_file_path, hist_names):
    """
    Load selected histograms (TH1/TH2/TH3/THnSparseD) from a ROOT file.

    Parameters
    ----------
    root_file_path : str
        Path to the ROOT file.
    hist_names : list or tuple of str
        Names of histograms to load from the file.

    Returns
    -------
    dict
        Dictionary {hist_name: histogram_object} with histograms detached from file.
    """
    if not os.path.isfile(root_file_path):
        raise FileNotFoundError(f"Cannot find ROOT file: {root_file_path}")

    # Open ROOT file
    f = ROOT.TFile.Open(root_file_path)
    if not f or f.IsZombie():
        raise IOError(f"Failed to open ROOT file: {root_file_path}")

    histograms = {}

    for name in hist_names:
        obj = f.Get(name)
        if not obj:
            print(f"Warning: Object '{name}' not found in file.")
            continue

        # Check for valid histogram-like type
        valid_types = (
            obj.InheritsFrom("TH1") or
            obj.InheritsFrom("TH2") or
            obj.InheritsFrom("TH3") or
            obj.InheritsFrom("THnSparse")
        )

        if not valid_types:
            print(f"Warning: Object '{name}' is not TH1/TH2/TH3/THnSparseD, skipping.")
            continue

        # Clone and detach from file
        h_clone = obj.Clone()
        h_clone.SetDirectory(0)

        histograms[name] = h_clone

    f.Close()
    return histograms

# Helper function to get input file path based on parameters
def input_file(category = "embedding", direction = "sum", eta_cut = 19, jet_selection = "jetId", systematics = None, trigger = "None"):

    #################### Parameters ###########################
    # - category: exp, embedding, pythia
    # - direction: pgoing, Pbgoing, sum
    # - eta_cut: integer value for eta cut in CM frame
    # - jet_selection: noCut, trkMax, jetId
    # - systematics:
    #    For experimental data:
    #       none, jeu_up, jeu_down, gplus, vtx1
    #    For MC:
    #       def, jerDef, jerUp, jerDown
    # - trigger: MB, Jet60, Jet80, Jet100, None - default, for MC
    ############################################################

    # Get user name
    username = getpass.getuser()
    # Base path
    base_path = Path(f"/Users/{username}/cernbox/ana/pPb8160")
    # Data directory
    data_dir = base_path / category

    # Handle special cases for direction
    if direction not in ("pgoing", "Pbgoing", "sum"):
        raise ValueError(f"Invalid direction '{direction}'. Must be 'pgoing', 'Pbgoing', or 'sum'.")

    # Default systematic choice for MC
    if category in ("embedding", "pythia") and systematics is None:
        systematics = "def"

    # File pattern logic
    if category == "embedding" or category == "pythia":
        # e.g. oEmbedding_pPb8160_def_ak4_jetId_eta19.root
        # or   Pbgoing/oEmbedding_Pbgoing_def_ak4_jetId_eta19.root

        category_name = "oEmbedding" if category == "embedding" else "oPythia"
        # For sum direction
        if direction == "sum":
            filename = f"{category_name}_pPb8160_{systematics}_ak4_{jet_selection}_eta{eta_cut}.root"
        else:
            # Direction-specific (pgoing or Pbgoing)
            filename = f"{direction}/{category_name}_{direction}_{systematics}_ak4_{jet_selection}_eta{eta_cut}.root"

        full_path = data_dir / filename
    else:
        # Experimental data
        # e.g. Jet100_pPb8160_ak4_jetId_eta19.root
        # or   Pbgoing/Jet100_Pbgoing_ak4_jetId_eta19.root

        if trigger is None or trigger == "None":
            raise ValueError("Trigger must be specified for experimental data.")

        if direction == "sum":
            if systematics is None:
                filename = f"{trigger}_pPb8160_ak4_{jet_selection}_eta{eta_cut}.root"
            else:
                filename = f"{trigger}_pPb8160_{systematics}_ak4_{jet_selection}_eta{eta_cut}.root"
        else:
            if systematics is None:
                filename = f"{direction}/{trigger}_{direction}_ak4_{jet_selection}_eta{eta_cut}.root"
            else:
                filename = f"{direction}/{trigger}_{direction}_{systematics}_ak4_{jet_selection}_eta{eta_cut}.root"
        full_path = data_dir / filename
    

    # Ensure filename is a Path and that the file exists before trying to open it
    input_file_name = Path(full_path)
    if not input_file_name.exists():
        raise FileNotFoundError(f"Input file not found: {input_file_name}")   
    return input_file_name

# -------------------------------------------------------
# Utility: project THnSparse with axis ranges
# -------------------------------------------------------
def project_sparse(sparse, axes=(0,1), ranges=None, name="proj"):
    """
    THnSparse projection that works in PyROOT.
    Supports:
      axes=(i)      → TH1D
      axes=(i,j)    → TH2D
      axes=(i,j,k)  → TH3D
    """

    ndim = sparse.GetNdimensions()

    # Save original axis ranges so we can restore them
    orig = [(sparse.GetAxis(i).GetFirst(),
             sparse.GetAxis(i).GetLast()) for i in range(ndim)]

    try:
        # Set cuts
        if ranges:
            for iax, (low, high) in ranges.items():
                ax = sparse.GetAxis(iax)
                first = ax.FindBin(low)
                last  = ax.FindBin(high*(1-1e-6))
                ax.SetRange(first,last)

        # ---------- ROOT-SAFE PROJECTION ----------
        if len(axes) == 1:
            proj = sparse.Projection(axes[0])

        elif len(axes) == 2:
            # order: Y, X  (ROOT convention)
            axY, axX = axes[0], axes[1]
            proj = sparse.Projection(axY, axX)

        elif len(axes) == 3:
            axZ, axY, axX = axes[0], axes[1], axes[2]
            proj = sparse.Projection(axZ, axY, axX)

        else:
            raise ValueError("Projection of >3 axes is not supported in this wrapper")

        # If ROOT misbehaves:
        if isinstance(proj, ROOT.THnSparse):
            raise RuntimeError(
                f"Projection returned THnSparse. Check axes={axes}."
            )

        proj.SetDirectory(0)
        proj.SetName(name.replace("-", "m"))
        return proj

    finally:
        # Restore original ranges
        for i,(first,last) in enumerate(orig):
            sparse.GetAxis(i).SetRange(first,last)

# -------------------------------------------------------
# FitSlicesY and return JES (mean) and JER (sigma) histos
# -------------------------------------------------------
def extract_jes_jer(h2):
    """
    h2 must be TH2D (y = pT_reco / pT_ref)
    Returns: (hMean, hSigma)
    """
    name = h2.GetName()
    h2.FitSlicesY()

    hMean  = ROOT.gDirectory.Get(f"{name}_1")
    hSigma = ROOT.gDirectory.Get(f"{name}_2")

    hMean.SetTitle("JES; ;JES")
    hSigma.SetTitle("JER; ;JER")

    return hMean, hSigma

#-------------------------------------------------------
# Format numeric range to safe name component
#-------------------------------------------------------
def fmt_range(low, high):
    """
    Convert numeric ranges to safe name components:
        -3.0 → m30
         0.8 → 08
        120  → 120
    """
    def f(x):
        s = f"{x}".replace('.', '_')
        if s.startswith('-'):
            s = 'm' + s[1:]
        return s
    return f(low), f(high)

#-------------------------------------------------------
# Get JES prefix from histogram name
#-------------------------------------------------------
def jes_prefix_from_name(hname):
    """Return clean short prefix ('hInclusive', 'hLead', 'hSubLead')."""
    if hname.startswith("hInclusive"):
        return "hInclusive"
    if hname.startswith("hLead"):
        return "hLead"
    if hname.startswith("hSubLead"):
        return "hSubLead"
    return "hUnknown"

#-------------------------------------------------------
# Draw TLatex slice label
#-------------------------------------------------------
def draw_slice_label(prefix, text, x=0.50, y=0.92):
    """
    Draw a TLatex label on the current pad.
    prefix: hInclusive/hLead/hSubLead → displayed as Inclusive/Lead/SubLead
    text: slice text, like "-1.6 < #eta < -0.8"
    """
    # Remove leading h if present
    clean_prefix = prefix[1:] if prefix.startswith("h") else prefix

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(22)
    latex.DrawLatex(x, y, f"{clean_prefix}: {text}")
