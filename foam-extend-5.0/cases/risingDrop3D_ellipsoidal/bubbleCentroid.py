"""Bubble centre of mass and rise velocity of an AMR case (run with pvpython).

The cell count changes at every adaptation, so a fixed-grid reader cannot be
used. ParaView's CellSize gives the exact cell volumes and CellCenters the
cell-centre coordinates, and the weighted sums are computed from the fetched
arrays:

    z_c = sum(alpha1 V z) / sum(alpha1 V),    v_c = dz_c/dt

with v_c by central differences. The last point is dropped (literature
convention for rising-bubble benchmarks).

Usage (normally called by postprocess.py):
    PYTHONNOUSERSITE=1 pvpython --force-offscreen-rendering \
        bubbleCentroid.py <case_dir> <out_csv>
"""

import glob
import os
import sys

import numpy as np


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
        OpenFOAMReader, CellSize, CellCenters, UpdatePipeline, servermanager,
    )

    case_dir, out_csv = sys.argv[1], sys.argv[2]
    reader = OpenFOAMReader(FileName=foam_file(case_dir))
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = ["alpha1"]
    reader.UpdatePipelineInformation()
    steps = [float(v) for v in (reader.TimestepValues or [])]
    print(f"{len(steps)} time steps in {case_dir}", flush=True)

    size = CellSize(Input=reader)
    size.ComputeLength = 0
    size.ComputeArea = 0
    size.ComputeVolume = 1
    size.ComputeVertexCount = 0

    centers = CellCenters(Input=size)
    centers.VertexCells = 0

    times, vols, zcs = [], [], []
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

        vol_sum = 0.0
        moment_sum = 0.0
        for block in blocks:
            # CellCenters outputs one point per input cell and carries the
            # original cell arrays as point data on those generated points.
            point_data = block.GetPointData()
            alpha_arr = point_data.GetArray("alpha1")
            volume_arr = point_data.GetArray("Volume")
            points = block.GetPoints()
            if alpha_arr is None or volume_arr is None or points is None:
                continue
            alpha = vtk_to_numpy(alpha_arr)
            volumes = vtk_to_numpy(volume_arr)
            xyz = vtk_to_numpy(points.GetData())
            weights = np.where(alpha > 1e-6, alpha * volumes, 0.0)
            vol_sum += float(weights.sum())
            moment_sum += float((weights * xyz[:, 2]).sum())

        if vol_sum <= 1e-12:
            print(f"  t = {t:.3f}: empty weighted volume, skipped",
                  flush=True)
            continue
        times.append(t)
        vols.append(vol_sum)
        zcs.append(moment_sum / vol_sum)

    t = np.array(times)
    vol = np.array(vols)
    zc = np.array(zcs)
    vc = np.zeros_like(zc)
    for i in range(1, len(t) - 1):
        dt = t[i + 1] - t[i - 1]
        if dt > 0:
            vc[i] = (zc[i + 1] - zc[i - 1]) / dt
    # Remove last point (literature convention).
    t, zc, vc, vol = t[:-1], zc[:-1], vc[:-1], vol[:-1]

    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    np.savetxt(out_csv, np.column_stack([t, zc, vc, vol]),
               header="t_s,z_c_m,v_c_m_s,volume_m3", delimiter=",",
               comments="")
    print(f"wrote {out_csv} ({len(t)} rows)", flush=True)


if __name__ == "__main__":
    main()
