/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dynamicT8codeFvMesh

Description
    Adaptive mesh refinement with t8code as the refinement backend: the
    t8code forest owns the mesh topology and the OpenFOAM polyMesh is rebuilt
    from the forest leaves after every adaptation.

    Each base cell becomes one t8code tree: a quad tree on quasi-2D meshes
    (any patch of type empty; the quad leaves are extruded back to one-layer
    hexahedra) and a hex tree on 3D meshes. Every refineInterval time steps:

    1. Cells are flagged where the point-interpolated refinement field lies in
       [lowerRefineLevel, upperRefineLevel), as in dynamicRefineFvMesh, and
       the band is protected from coarsening together with nBufferLayers
       cell layers around it and, optionally, every cell whose centre lies
       within bufferDistance of a band cell. Only the band itself is
       refined.
    2. t8code refines the band cells by one level, up to maxRefinement, and
       coarsens the families with no protected member whose smallest value
       is below unrefineLevel (the native rule; unrefineLevel above the
       field range, e.g. 10 for alpha1, coarsens both phases).
    3. alpha1, U and pd are interpolated in t8code callbacks. The polyMesh is
       rebuilt through directTopoChange, either from scratch or, with
       incrementalRebuild, only where t8code changed the forest. The
       resulting mapPolyMesh maps all remaining fields, and the fluxes in
       correctFluxes are recomputed only on faces without an old master.

    Each call of update() logs one line with the cost of the whole mesh
    update: `AMRmeshUpdate: CPU <s> s  Clock <s> s  changed <0|1>`.

    Conversion between the two representations:
    - OF -> t8: initForest(), copyFieldsToForest()
    - t8 -> OF: forestToOpenFOAM(), forestToOpenFOAMIncremental(),
      copyFieldsFromForest()

    Lifetime rule (avoids the dangling user_data segfault): the per-element
    field arrays are std::vector members of this class, so they outlive every
    adapt call; user_data is cleared before any t8_forest_unref.

    Limitations: serial runs only; the interpolated fields are fixed to the
    VOF set (the refinement field, U and pd); maxCells is required for
    compatibility with dynamicRefineFvMesh but not used yet.

Usage
    Load the library in system/controlDict:
    \verbatim
    libs
    (
        "libt8codeDynamicFvMesh.so"
    );
    \endverbatim

    and select the mesh in constant/dynamicMeshDict:
    \verbatim
    dynamicFvMesh   dynamicT8codeFvMesh;

    dynamicT8codeFvMeshCoeffs
    {
        field               alpha1;     // refinement field
        refineInterval      5;          // adapt every N steps (<= 0: never)
        lowerRefineLevel    0.001;      // refine where field is in
        upperRefineLevel    0.999;      //   [lower, upper)
        unrefineLevel       10;         // coarsen where min field < this
        maxRefinement       2;          // maximum refinement level
        maxCells            2000000;    // required, not used yet
        nBufferLayers       1;          // optional (default 1)
        bufferDistance      0;          // optional (default 0: off)
        incrementalRebuild  true;       // optional (default false)
        correctFluxes       ((phi U));  // fluxes recomputed after mapping
    }
    \endverbatim

    Testing only: with the environment variable T8CHECK_MAP set, the
    leaf-to-cell map is cross-checked against a full findCell search after
    every rebuild.

SourceFiles
    dynamicT8codeFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicT8codeFvMesh_H
#define dynamicT8codeFvMesh_H

#include "dynamicFvMesh.H"

#include <vector>

namespace Foam
{

// Forward declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                         Struct T8FieldData Declaration
\*---------------------------------------------------------------------------*/

//- Plain (t8code-free) view of the per-element fields handed to the C
//  callbacks
struct T8FieldData
{
    //- Refinement field (alpha1), one value per element
    double* alpha1;

    //- Velocity, 3 components per element, interleaved
    double* U;

    //- Dynamic pressure, one value per element
    double* pd;

    //- Number of elements in the arrays
    int     nElements;
};


/*---------------------------------------------------------------------------*\
                     Class dynamicT8codeFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicT8codeFvMesh
:
    public dynamicFvMesh
{
    // Private data

        // Run-time controls (read from dynamicT8codeFvMeshCoeffs)

            //- Adapt every refineInterval time steps (<= 0: never)
            label refineInterval_;

            //- Name of the refinement field
            word  fieldName_;

            //- Lower bound of the refinement band
            scalar lowerRefineLevel_;

            //- Upper bound of the refinement band
            scalar upperRefineLevel_;

            //- A family is coarsened only if its smallest value is below
            //  this level (dynamicRefineFvMesh semantics)
            scalar unrefineLevel_;

            //- Maximum refinement level of a base cell
            label maxRefinement_;

            //- Read for compatibility with dynamicRefineFvMesh; not used
            label maxCells_;

        //- Number of cell layers around the interface band protected from
        //  coarsening (as in dynamicRefineFvMesh). Read from
        //  dynamicT8codeFvMeshCoeffs::nBufferLayers (default 1).
        label nBufferLayers_;

        //- Cells whose centre lies within this distance of a band cell
        //  centre are also protected from coarsening (0: off). Read from
        //  dynamicT8codeFvMeshCoeffs::bufferDistance (default 0).
        scalar bufferDistance_;

        //- Fluxes recomputed after mapping, as (flux velocity) pairs
        List<Pair<word>> correctFluxes_;

        //- Per-leaf refine flag (interface band only); lifetime spans the
        //  adapt call that reads it
        std::vector<double> refineFlag_;

        //- Per-leaf protect flag (interface band dilated by nBufferLayers_):
        //  these cells are not coarsened, as in dynamicRefineFvMesh
        std::vector<double> protectFlag_;

        //- Per NEW-cell -> OLD-cell map, built from the t8 parent<->child
        //  relationship (t8_forest_iterate_replace). Used as the directTopoChange
        //  masterCell for field mapping, replacing an expensive geometric findCell.
        std::vector<label> masterCell_;

        //- Incremental rebuild controls (dirty-region only, like hexRef8).
        //  When true, only the cells t8 changed are reconstructed; the rest of
        //  the OF mesh is retained. From-scratch rebuild is the fallback.
        bool incrementalRebuild_;

        //- forest leaf -> OF cell index, maintained across rebuilds (identity
        //  after a from-scratch rebuild; updated through the mapPolyMesh after
        //  an incremental one). Size = current nCells / nLeaf.
        std::vector<label> leafToCell_;

        //- per NEW leaf: 1 if it is part of a refine/coarsen (dirty), else 0.
        //  Built by the iterate_replace callback alongside masterCell_.
        std::vector<unsigned char> leafDirty_;

        //- per OLD leaf: 1 if removed by a refine/coarsen (its OF cell must be
        //  deleted in the incremental rebuild), else 0.
        std::vector<unsigned char> oldDirty_;

        //- Time bookkeeping
        label curTimeIndex_;

        // Geometry of the collapsed (empty) direction for the 2D extrusion

            //- True when the mesh has an empty patch (quad trees)
            bool planar2D_;

            //- Collapsed direction of a quasi-2D mesh
            direction emptyDir_;

            //- Extent of the mesh along the collapsed direction
            scalar zMin_;
            scalar zMax_;

        //- Number of base (level-0) trees == number of base cells
        label nBaseTrees_;

        // t8code objects (opaque; real types only in the .C)

            //- Current forest (t8_forest_t)
            void* forest_;

            //- Coarse mesh (t8_cmesh_t), owned by the forest once set
            void* cmesh_;

            //- Refinement scheme (const t8_scheme_c*)
            void* scheme_;

            //- True after initForest() has built the level-0 forest
            bool  forestInitialized_;

        // Per-element field storage (lifetime spans every adapt call)

            //- Refinement field per element
            std::vector<double> alpha1Data_;

            //- Velocity per element (3 interleaved components)
            std::vector<double> uData_;

            //- Dynamic pressure per element
            std::vector<double> pdData_;

        //- Field view passed to the C callbacks (points into the vectors above)
        T8FieldData fieldData_;


    // Private Member Functions

        //- Read dynamicT8codeFvMeshCoeffs from constant/dynamicMeshDict
        void readDict();

        //- OF -> t8 : build the cmesh (quad trees in 2D, hex trees in 3D)
        //  and the level-0 forest from this mesh
        void initForest();

        //- OF -> t8 : copy alpha1/U/pd into the per-element vectors + user_data
        void copyFieldsToForest();

        //- Drive one t8code adapt (decision) + balance + field interpolation
        bool runT8codeAdapt();

        //- t8 -> OF : rebuild the WHOLE polyMesh from the forest leaves
        bool forestToOpenFOAM();

        //- t8 -> OF : rebuild only the cells t8 changed (+ finer-neighbour halo)
        bool forestToOpenFOAMIncremental();

        //- Re-derive leafToCell_ (forest leaf <-> OF cell) by sorting the
        //  leaf centroids and the cell centres in the same order; pairs
        //  that fail the distance check fall back to findCell
        void rebuildLeafToCell();

        //- t8 -> OF : write the interpolated per-element fields back
        void copyFieldsFromForest();

        //- Body of update() (forest init, adapt and OF rebuild), without
        //  the once-per-time-step guard and the AMRmeshUpdate log line
        bool updateTopology();

        //- Recompute the listed fluxes from the mapped velocity, but only
        //  on faces without an old master. Faces with a valid master keep
        //  the mapped (MULES-limited) flux, exactly like
        //  dynamicRefineFvMesh. Wholesale overwrite would discard the
        //  limiter at every adapt and feed dispersive interface errors.
        void correctMappedFluxes(const mapPolyMesh& map);

        //- Disallow default bitwise copy
        dynamicT8codeFvMesh(const dynamicT8codeFvMesh&);
        void operator=(const dynamicT8codeFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicT8codeFvMesh");


    // Constructors

        //- Construct from IOobject (reads dynamicT8codeFvMeshCoeffs; the
        //  forest is built on the first call of update())
        explicit dynamicT8codeFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicT8codeFvMesh();


    // Member Functions

        //- Adapt the mesh (once per time step, every refineInterval steps)
        //  and log the AMRmeshUpdate line. Returns true if the topology
        //  changed.
        virtual bool update();
};


} // End namespace Foam

#endif

// ************************************************************************* //
