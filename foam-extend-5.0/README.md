# foamT8codeAMR — foam-extend-5.0

Adaptive mesh refinement (AMR) for
[foam-extend-5.0](https://sourceforge.net/projects/foam-extend/) using
[t8code](https://github.com/DLR-AMR/t8code) as the refinement backend,
replacing the native `hexRef8` engine.

## Contents

|Path|Description|
|---|---|
|`src/t8codeDynamicFvMesh/`|`dynamicT8codeFvMesh` library (t8code forest owns the topology)|
|`applications/solvers/interT8codeDyMFoam/`|interDyMFoam copy linked against the t8code library|
|`cases/damBreakValidation_AMR_t8code/`|2D dam-break validation case with experimental gauge data|
|`cases/risingDrop_ellipsoidal_AMR_3s_t8code/`|Ellipsoidal rising drop case with reference data|
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

## Usage

Select the mesh in `constant/dynamicMeshDict`:

```cpp
dynamicFvMesh   dynamicT8codeFvMesh;
```

and see the example cases for the full coefficients dictionary. Each
case is run with `./Allrun` and cleaned with `./Allclean`. The `Allrun`
scripts source a local foam-extend path; adjust it to your installation.
