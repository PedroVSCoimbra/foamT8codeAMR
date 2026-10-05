# foamT8codeAMR — foam-extend-5.0

Adaptive mesh refinement (AMR) for
[foam-extend-5.0](https://sourceforge.net/projects/foam-extend/) using
[t8code](https://github.com/DLR-AMR/t8code) as the refinement backend,
replacing the native `hexRef8` engine.

## Contents

|Path|Description|
|---|---|
|`src/t8codeDynamicFvMesh/`|`dynamicT8codeFvMesh` library (t8code forest owns the topology)|
|`applications/solvers/interT8codeDyMFoam/`|interDyMFoam built as a separate executable linked against the t8code library|
|`cases/damBreak2D/`|2D dam break, gauge heights vs. Zhainakov & Kurbanaliev (2013)|
|`cases/risingDrop3D_ellipsoidal/`|3D rising drop (ellipsoidal regime) vs. Silva et al. (2023)|
|`docs/`|Architecture notes (in Portuguese)|

The library and solver are built out of tree: no file of foam-extend-5.0
itself is modified.

## Requirements

- foam-extend-5.0 (C++14 toolchain, tested with the gcc-11 port)
- t8code built with MPI (tested with `v4.0.0-26.06-79-g43370a122`),
  by default expected at `$HOME/t8code` with its CMake build in
  `$HOME/t8code/build`

## Build

```bash
source /path/to/foam-extend-5.0/etc/bashrc
./Allwmake                              # t8code in $HOME/t8code
T8CODE_DIR=/path/to/t8code ./Allwmake   # custom t8code location
```

The library goes to `$FOAM_LIBBIN` and the solver to `$FOAM_APPBIN`.
The t8code, sc and p4est library directories are embedded as RPATH.

`interT8codeDyMFoam` contains only `interT8codeDyMFoam.C`, a copy of the
upstream `interDyMFoam.C`. Its other headers are read from
`$FOAM_SOLVERS/multiphase/interDyMFoam`, so the solver algorithm is
exactly that of `interDyMFoam`. It exists as a separate executable only to
link `libt8codeDynamicFvMesh` without touching the original solver.

## Usage

Select the mesh in `constant/dynamicMeshDict`:

```cpp
dynamicFvMesh   dynamicT8codeFvMesh;
```

and see the example cases for the full coefficients dictionary. Each
case has its own README with the setup, the reference data and how they
were obtained. Cases are run with `./Allrun`, post-processed with
`python3 postprocess.py` and cleaned with `./Allclean`. The `Allrun`
scripts source a local foam-extend path; adjust it to your installation.

## Log output

At every time step the mesh prints one line with the cost of the whole
mesh update, adapting or not:

```text
AMRmeshUpdate: CPU <s> s  Clock <s> s  changed <0|1>
```

On adaptation steps, a breakdown is printed as well
(`dynamicT8codeFvMesh::update timing [ms]: t8 adapt = ..., OF rebuild = ...`).
