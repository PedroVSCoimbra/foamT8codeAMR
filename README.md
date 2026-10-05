# foamT8codeAMR

Adaptive mesh refinement (AMR) for OpenFOAM-family codes using
[t8code](https://github.com/DLR-AMR/t8code) as the refinement backend.

Each supported OpenFOAM fork/version lives in its own directory, with its
own library, solvers, example cases and build scripts:

|Directory|Target|Status|
|---|---|---|
|[`foam-extend-5.0/`](foam-extend-5.0/README.md)|foam-extend-5.0|In development|

All code is built out of tree: no file of the host OpenFOAM installation
is modified.

## License

GNU General Public License v3.0 or later (see `COPYING`).
