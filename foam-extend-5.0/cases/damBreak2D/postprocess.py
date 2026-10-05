#!/usr/bin/env python3
"""Water levels at the gauges x1 = 2.725 m and x2 = 2.228 m.

Compares the free-surface height h(t) of the AMR run with the markers of
Zhainakov & Kurbanaliev (2013), Fig. 4 (see README.md). The present-work
curve is computed by gaugeHeights.py (pvpython) and cached in
postProcessing/gauge_heights.csv; delete that file to recompute.

Usage (from the case directory, after ./Allrun):
    python3 postprocess.py
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

CASE_DIR = Path(__file__).resolve().parent
OUT_DIR = CASE_DIR / "postProcessing"
REF_DIR = CASE_DIR / "reference"
PRESENT_CSV = OUT_DIR / "gauge_heights.csv"

# (reference file, column in PRESENT_CSV)
GAUGES = [
    ("zhainakov2013_fig4a_x2.725m.csv", 1),
    ("zhainakov2013_fig4b_x2.228m.csv", 2),
]

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "text.usetex": True,
})


def load_reference(path):
    t, h = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["t_s"]))
            h.append(float(row["h_m"]))
    return np.array(t), np.array(h)


def present_curve():
    """Rows (t, h1, h2) of the present run, computing them if not cached."""
    if not PRESENT_CSV.exists():
        script = CASE_DIR / "gaugeHeights.py"
        env = dict(os.environ, PYTHONNOUSERSITE="1")
        subprocess.run(
            ["pvpython", "--force-offscreen-rendering", str(script),
             str(CASE_DIR), str(PRESENT_CSV)],
            check=True, env=env,
        )
    return np.loadtxt(PRESENT_CSV, skiprows=1, delimiter=",", ndmin=2)


def main():
    try:
        present = present_curve()
    except (OSError, subprocess.CalledProcessError) as err:
        sys.exit(f"Could not compute the gauge heights: {err}")

    fig, axes = plt.subplots(2, 1, figsize=(16 / 2.54, 12 / 2.54),
                             sharex=True)
    for ax, (ref_file, col) in zip(axes, GAUGES):
        ref_t, ref_h = load_reference(REF_DIR / ref_file)
        # Same time window as the reference; last point dropped
        # (literature convention)
        keep = present[:, 0] <= ref_t[-1]
        t, h = present[keep, 0][:-1], present[keep, col][:-1]
        ax.plot(ref_t, ref_h, "o", markerfacecolor="none",
                markeredgecolor="k", markersize=6,
                label=r"Zhainakov \& Kurbanaliev (2013)")
        ax.plot(t, h, "-", lw=1.5, label="Present (foam-extend + t8code)")
        ax.set_ylabel(r"$h$ [m]")
        ax.grid(True, alpha=0.3)
    axes[0].legend(loc="upper left")
    axes[1].set_xlabel(r"$t$ [s]")
    fig.tight_layout()

    OUT_DIR.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"water_levels.{ext}", dpi=200)
    print(f"Saved {OUT_DIR / 'water_levels.pdf'} (+ .png)")


if __name__ == "__main__":
    main()
