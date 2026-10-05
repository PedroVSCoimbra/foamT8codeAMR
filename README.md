# foamT8codeAMR

Adaptive mesh refinement (AMR) for OpenFOAM-family codes using
[t8code](https://github.com/DLR-AMR/t8code) as the refinement backend.

The code is organized in one directory per OpenFOAM fork/version (e.g.
`foam-extend-5.0/`, `openfoam-9/`, `openfoam-v26xx/`), each with its own
library, solvers, example cases, build scripts and README. Currently only
[`foam-extend-5.0/`](foam-extend-5.0/README.md) is available.

All code is built out of tree: no file of the host OpenFOAM installation
is modified.

## License

GNU General Public License v3.0 or later (see `COPYING`).
