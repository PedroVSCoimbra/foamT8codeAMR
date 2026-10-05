"""Free-surface height at the gauges of an AMR case (run with pvpython).

The mesh changes at every adaptation, so the cell centres are read for each
time step by ParaView (CellCenters) instead of parsing polyMesh in Python.
At each gauge x_g the height is the highest cell centre with alpha1 > 0.5
inside the strip |x - x_g| < 0.05 m (0 when the strip is dry).

Usage (normally called by postprocess.py):
    PYTHONNOUSERSITE=1 pvpython --force-offscreen-rendering \
        gaugeHeights.py <case_dir> <out_csv>
"""

import glob
import os
import sys

import numpy as np

GAUGES = (2.725, 2.228)  # x1, x2 [m]
HALF_WIDTH = 0.05        # strip half-width [m]


def foam_file(case_dir):
    existing = glob.glob(os.path.join(case_dir, "*.foam"))
    if existing:
        return existing[0]
    foam = os.path.join(case_dir, "case.foam")
    open(foam, "a").close()
    return foam


def main():
    # pylint: disable=import-outside-toplevel
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from paraview.simple import (
        OpenFOAMReader, CellCenters, UpdatePipeline, servermanager,
    )

    case_dir, out_csv = sys.argv[1], sys.argv[2]
    reader = OpenFOAMReader(FileName=foam_file(case_dir))
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = ["alpha1"]
    reader.UpdatePipelineInformation()
    steps = [float(v) for v in (reader.TimestepValues or [])]
    print(f"{len(steps)} time steps in {case_dir}", flush=True)

    centers = CellCenters(Input=reader)
    centers.VertexCells = 0

    rows = []
    for t in steps:
        UpdatePipeline(time=t, proxy=centers)
        data = servermanager.Fetch(centers)
        blocks = []
        if hasattr(data, "GetNumberOfBlocks"):
            for block_i in range(data.GetNumberOfBlocks()):
                block = data.GetBlock(block_i)
                if block is not None:
                    blocks.append(block)
        else:
            blocks.append(data)

        xyz_list, alpha_list = [], []
        for block in blocks:
            alpha_arr = block.GetPointData().GetArray("alpha1")
            points = block.GetPoints()
            if alpha_arr is None or points is None:
                continue
            xyz_list.append(vtk_to_numpy(points.GetData()))
            alpha_list.append(vtk_to_numpy(alpha_arr))
        if not xyz_list:
            continue
        xyz = np.concatenate(xyz_list)
        alpha = np.concatenate(alpha_list)

        heights = []
        for xg in GAUGES:
            wet = (np.abs(xyz[:, 0] - xg) < HALF_WIDTH) & (alpha > 0.5)
            heights.append(float(xyz[wet, 1].max()) if wet.any() else 0.0)
        rows.append([t] + heights)

    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    np.savetxt(out_csv, np.array(rows), header="t_s,h1_m,h2_m",
               delimiter=",", comments="")
    print(f"wrote {out_csv} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
