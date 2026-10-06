# foamT8codeAMR

Adaptive mesh refinement (AMR) for OpenFOAM-family codes using
[t8code](https://github.com/DLR-AMR/t8code) (spoken as "tetcode") as the
refinement backend. The AMR is a run-time selectable `dynamicFvMesh`
library, loaded through `libs` in `system/controlDict` and used with the
standard dynamic-mesh solvers.

## Directory organization

```text
foamT8codeAMR/
├── doc/        architecture notes, Doxygen configuration
├── fe50/       foam-extend-5.0
│   ├── src/        Allwmake, libs/
│   └── tutorials/  validation cases with reference data
├── of90/       OpenFOAM 9 (not available yet)
└── ofv2606/    OpenFOAM v2606 (not available yet)
```

Each version directory is self-contained and has its own README. Only
[`fe50/`](fe50/README.md) is available for now; the other version
directories are placeholders. After cloning, the directories of the
versions you do not use can be removed.

All code is built out of tree: no file of the host OpenFOAM installation
is modified.

## Documentation

The source documentation is generated with Doxygen (requires `doxygen`,
`python3` and, for the class graphs, `graphviz`):

```bash
doc/Doxygen/Allwmake    # then open doc/Doxygen/html/index.html
```

## License

GNU General Public License v3.0 or later (see `COPYING`).
