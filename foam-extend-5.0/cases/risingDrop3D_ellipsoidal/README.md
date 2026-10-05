# risingDrop3D_ellipsoidal

3D buoyancy-driven rising drop, test case 1 (ellipsoidal regime), run with
`interT8codeDyMFoam` and t8code AMR. The bubble centre of mass and rise
velocity are compared with the solvers reported by Silva et al. (2023).

## Setup

|Item|Value|
|---|---|
|Domain (x × y × z)|1 × 1 × 2, walls on all sides|
|Drop|sphere, radius 0.25, centre (0.5, 0.5, 0.5)|
|Outer fluid (`phase2`)|ρ₁ = 1000, μ₁ = 10|
|Drop (`phase1`)|ρ₂ = 100, μ₂ = 1|
|Surface tension|σ = 24.5|
|Gravity|(0, 0, −0.98)|
|Dimensionless|Re = 35, Bo = 10|
|Base mesh|25 × 25 × 50|
|AMR|`maxRefinement 2` (Δx/4 at the interface), `refineInterval 5`, incremental rebuild|
|Time|t = 0 to 3, `maxCo 0.5`, written every 0.05|

Physical parameters follow Table 3 (case 1) of the reference. Quantities
are dimensionless, as in the benchmark. In foam-extend, `nu` is set to
μ/ρ = 0.01 for both phases.

## Running

```bash
./Allrun              # blockMesh, setFields, interT8codeDyMFoam
python3 postprocess.py
./Allclean
```

`postprocess.py` calls `bubbleCentroid.py` through `pvpython` (ParaView)
to compute, at each written time,

- the centre of mass z_c = Σ α V z / Σ α V, with exact AMR cell volumes,
- the rise velocity v_c = dz_c/dt, by central differences.

It caches the result in `postProcessing/bubble_data.csv` and writes
`postProcessing/bubble_evolution.pdf` (+ `.png`). Following the
literature convention, the last time point is dropped. Delete the CSV to
recompute.

## Reference data

`reference/silva2023_fig9_case1.csv`

|Column|Content|
|---|---|
|`t`|time [-]|
|`z_<solver>`|bubble centre of mass z_c [-] (Fig. 9a)|
|`w_<solver>`|bubble rise velocity w [-] (Fig. 9b)|

Solvers: `HiG-Flow`, `OpenFOAM`, `DROPS`, `FEEL`, `NaSt3D`. Empty cells
are points that could not be read from the figure.

- **Source:** A. T. G. da Silva, C. Fernandes, J. Organista, L. Souza and
  A. Castelo, "Validation of HiG-Flow Software for Simulating Two-Phase
  Flows with a 3D Geometric Volume of Fluid Algorithm", *Mathematics*
  11(18), 3900, 2023. <https://doi.org/10.3390/math11183900>. Figure 9:
  "3D buoyancy-driven rising drop: Test case 1 (ellipsoidal regime)".
- **Origin of each curve:**
  - HiG-Flow: Silva et al. (2023).
  - FEEL: Metivet et al., *Mathematics* 6, 203, 2018.
  - OpenFOAM, DROPS and NaSt3D: Adelsberger et al., "3D incompressible
    two-phase flow benchmark computations for rising droplets", WCCM XI,
    Barcelona, 2014.

  All five are as plotted in Silva et al. (2023), Fig. 9.
- **Accuracy:** values were read from the published figure by manual
  digitization and are approximate.
