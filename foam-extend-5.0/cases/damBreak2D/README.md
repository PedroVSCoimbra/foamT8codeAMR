# damBreak2D

2D dam break in a closed tank, run with `interT8codeDyMFoam` and t8code
AMR. The free-surface height at two gauges is compared with the data
reported by Zhainakov & Kurbanaliev (2013).

## Setup

|Item|Value|
|---|---|
|Tank (x × y)|3.22 m × 2.0 m, 2D (one cell in z, `frontAndBack` empty)|
|Water column|1.2 m × 0.6 m in the lower-left corner, at rest at t = 0|
|Fluids|water ρ = 998.2 kg/m³, ν = 1e-6 m²/s; air ρ = 1.225 kg/m³, ν = 1.48e-5 m²/s|
|Gravity|(0, −9.81, 0) m/s²|
|Turbulence|RAS k-ε|
|Gauges|x₁ = 2.725 m and x₂ = 2.228 m|
|Base mesh|160 × 100 × 1|
|AMR|`maxRefinement 2` (Δx/4 at the interface), `refineInterval 5`|
|Time|t = 0 to 2.5 s, `maxCo 0.1`, written every 0.01 s|

The geometry follows Section 3 of the reference: a tank 3.22 m long and
2 m high with a 1.2 m × 0.6 m water column, modelled in 2D as in the
paper. The paper's pressure probe P on the right wall is not used here.

## Running

```bash
./Allrun              # blockMesh, setFields, interT8codeDyMFoam
python3 postprocess.py
./Allclean
```

`postprocess.py` calls `gaugeHeights.py` through `pvpython` (ParaView),
which reads the adaptive mesh at each written time. The height h(t) at
each gauge is the highest cell centre with `alpha1 > 0.5` within ±5 cm
of the gauge. The result is cached in `postProcessing/gauge_heights.csv`
(delete it to recompute), and the figure is written to
`postProcessing/water_levels.pdf` (+ `.png`).

## Reference data

|File|Content|
|---|---|
|`reference/zhainakov2013_fig4a_x2.725m.csv`|h(t) at x₁ = 2.725 m (Fig. 4a)|
|`reference/zhainakov2013_fig4b_x2.228m.csv`|h(t) at x₂ = 2.228 m (Fig. 4b)|

Columns: `t_s` (time [s]), `h_m` (free-surface height above the tank
floor [m]).

- **Source:** A. Zh. Zhainakov and A. Y. Kurbanaliev, "Verification of
  the open package OpenFOAM on dam break problems", *Thermophysics and
  Aeromechanics* 20(4), 434–461, 2013. Figure 4: "Water column height in
  different sections. x₁ = 2.725 m (a) and x₂ = 2.228 m (b)".
- **What was digitized:** the markers of Fig. 4, which the paper
  describes as experimental data. The paper takes the test configuration
  from Zhou, De Kat & Buchner (1999), its reference [5]. The solid lines
  (the paper's OpenFOAM 1.7.1 results) were not digitized.
- **Accuracy:** values were read from the published figure by manual
  digitization; expect errors of the order of the marker size (~1 cm in
  h, ~0.02 s in t).
- **Correction:** the first markers of both panels lie on h = 0 (dry bed
  before the front arrives). Digitization gave them small negative values
  (down to −0.009 m), which were set to 0 (6 points in Fig. 4a, 3 in
  Fig. 4b). No other value was modified.
- The paper notes that agreement at x₁ is satisfactory only up to
  t ≈ 1.5 s, when the reflected wave reaches the gauge.
