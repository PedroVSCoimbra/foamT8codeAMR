# foamT8codeAMR — foam-extend-5.0 (`fe50`)

Adaptive mesh refinement (AMR) for
[foam-extend-5.0](https://sourceforge.net/projects/foam-extend/) using
[t8code](https://github.com/DLR-AMR/t8code) (spoken as "tetcode") as the
refinement backend, replacing the native `hexRef8` engine.

The AMR is a run-time selectable `dynamicFvMesh` in a shared library: it
works with the standard solvers (e.g. `interDyMFoam`), with no custom
solver and no change to foam-extend-5.0 itself.

## Contents

|Path|Description|
|---|---|
|`src/libs/t8codeDynamicFvMesh/`|`dynamicT8codeFvMesh` library (t8code forest owns the topology)|
|`src/Allwmake`, `src/Allwclean`|Build / clean scripts|
|`tutorials/damBreak2D/`|2D dam break, gauge heights vs. Zhainakov & Kurbanaliev (2013)|
|`tutorials/risingDrop3D_ellipsoidal/`|3D rising drop (ellipsoidal regime) vs. Silva et al. (2023)|

Architecture notes are in the top-level [`doc/`](../doc/) directory.

## Requirements

- foam-extend-5.0 (C++14 toolchain, tested with the gcc-11 port)
- t8code built with MPI (tested with `v4.0.0-26.06-79-g43370a122`),
  by default expected at `$HOME/t8code` with its CMake build in
  `$HOME/t8code/build`

## Build

```bash
source /path/to/foam-extend-5.0/etc/bashrc
cd src
./Allwmake                              # t8code in $HOME/t8code
T8CODE_DIR=/path/to/t8code ./Allwmake   # custom t8code location
```

The library goes to `$FOAM_LIBBIN`. The t8code, sc and p4est library
directories are embedded as RPATH, so no `LD_LIBRARY_PATH` change is
needed.

## Usage

Load the library in `system/controlDict`:

```cpp
libs
(
    "libt8codeDynamicFvMesh.so"
);
```

select the mesh in `constant/dynamicMeshDict`:

```cpp
dynamicFvMesh   dynamicT8codeFvMesh;
```

and run a dynamic-mesh solver such as `interDyMFoam`. See the tutorials
for the full `dynamicT8codeFvMeshCoeffs` dictionary. Each tutorial has its
own README with the setup, the reference data and how they were obtained.
Tutorials are run with `./Allrun`, post-processed with
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
