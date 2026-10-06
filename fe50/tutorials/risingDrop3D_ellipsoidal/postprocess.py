#!/usr/bin/env python3
"""Rising drop, test case 1 (ellipsoidal regime): z_c(t) and v_c(t).

Compares the bubble centre of mass and rise velocity of the AMR run with
the curves digitized from Silva et al. (2023), Fig. 9 (see README.md).
The present-work curve is computed by bubbleCentroid.py (pvpython, exact
AMR cell volumes) and cached in postProcessing/bubble_data.csv; delete
that file to recompute.

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
REF_CSV = CASE_DIR / "reference" / "silva2023_fig9_case1.csv"
PRESENT_CSV = OUT_DIR / "bubble_data.csv"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "text.usetex": True,
})

# Reference solvers: (label, z column, w column, colour)
REFERENCES = [
    ("HiG-Flow", "z_HiG-Flow", "w_HiG-Flow", "#1f77b4"),
    ("OpenFOAM", "z_OpenFOAM", "w_OpenFOAM", "#ff7f0e"),
    ("DROPS", "z_DROPS", "w_DROPS", "#2ca02c"),
    ("FEEL", "z_FEEL", "w_FEEL", "#d62728"),
    ("NaSt3D", "z_NaSt3D", "w_NaSt3D", "#9467bd"),
]


def load_reference(path):
    """Reference columns keyed by name (NaN where the digitization has gaps)."""
    t, data = [], {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t.append(float(row["t"]))
            for key, val in row.items():
                if key != "t":
                    val = val.strip()
                    data.setdefault(key, []).append(
                        float(val) if val else np.nan)
    return np.array(t), {k: np.array(v) for k, v in data.items()}


def present_curve():
    """(t, z_c, v_c) of the present run, computing it if not cached."""
    if not PRESENT_CSV.exists():
        script = CASE_DIR / "bubbleCentroid.py"
        env = dict(os.environ, PYTHONNOUSERSITE="1")
        subprocess.run(
            ["pvpython", "--force-offscreen-rendering", str(script),
             str(CASE_DIR), str(PRESENT_CSV)],
            check=True, env=env,
        )
    d = np.loadtxt(PRESENT_CSV, skiprows=1, delimiter=",", ndmin=2)
    return d[:, 0], d[:, 1], d[:, 2]


def main():
    ref_t, ref = load_reference(REF_CSV)
    try:
        t, zc, vc = present_curve()
    except (OSError, subprocess.CalledProcessError) as err:
        sys.exit(f"Could not compute the bubble centroid: {err}")

    ip = int(np.argmax(vc))
    print(f"Peak rise velocity: {vc[ip]:.3f} m/s at t = {t[ip]:.2f} s")
    print(f"Final centre of mass: {zc[-1]:.3f} m at t = {t[-1]:.2f} s")

    fig, (az, av) = plt.subplots(1, 2, figsize=(16 / 2.54, 3.6),
                                 sharex=True)
    for label, zkey, wkey, colour in REFERENCES:
        az.plot(ref_t, ref[zkey], "-", color=colour, lw=1.2)
        av.plot(ref_t, ref[wkey], "-", color=colour, lw=1.2)
    az.plot(t, zc, "k--", lw=1.8)
    av.plot(t, vc, "k--", lw=1.8)

    for ax in (az, av):
        ax.set_xlabel(r"$t$ [s]")
        ax.set_xlim(0, 3)
        ax.grid(True, alpha=0.3)
    az.set_ylabel(r"$z_c$ [m]")
    av.set_ylabel(r"$v_c$ [m/s]")
    az.set_ylim(0.4, 1.5)
    av.set_ylim(0, 0.42)

    handles = [plt.Line2D([], [], color=c, lw=1.2, label=lbl)
               for lbl, _, _, c in REFERENCES]
    handles.append(plt.Line2D([], [], color="k", ls="--", lw=1.8,
                              label="Present (foam-extend + t8code)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=True,
               bbox_to_anchor=(0.5, 0.015), borderaxespad=0.0)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.96, bottom=0.36,
                        wspace=0.3)

    OUT_DIR.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"bubble_evolution.{ext}", dpi=200)
    print(f"Saved {OUT_DIR / 'bubble_evolution.pdf'} (+ .png)")


if __name__ == "__main__":
    main()
