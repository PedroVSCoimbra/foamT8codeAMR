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
    A dynamicFvMesh in which t8code OWNS the mesh topology.

    Phase 1 (serial, 2D quasi-mesh): the in-plane footprint of every base cell
    becomes one quad tree of a t8code cmesh; the forest decides refine/coarsen
    from the alpha1 field; the OpenFOAM polyMesh is then REBUILT from the forest
    leaves (each quad extruded back to a one-layer hex) and the interpolated
    fields are written back.

    Conversion systems:
      - OF -> t8 : initForest()  (cmesh = 1 quad tree per base cell)
                   copyFieldsToForest()
      - t8 -> OF : forestToOpenFOAM()  (leaves -> points/faces/cells/patches)
                   copyFieldsFromForest()

    Lifetime rule (avoids the dangling-user_data segfault): the per-element
    field arrays are std::vector members of this class, so they outlive every
    adapt call; user_data is cleared before any t8_forest_unref.

\*---------------------------------------------------------------------------*/

#ifndef dynamicT8codeFvMesh_H
#define dynamicT8codeFvMesh_H

#include "dynamicFvMesh.H"

#include <vector>

namespace Foam
{

// Forward declarations
class mapPolyMesh;

// Plain (t8code-free) view of the per-element fields handed to the C callbacks
struct T8FieldData
{
    double* alpha1;
    double* U;          // 3 components per element, interleaved
    double* pd;
    int     nElements;
};

class dynamicT8codeFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Run-time controls (read from dynamicT8codeFvMeshCoeffs)
        label refineInterval_;
        word  fieldName_;
        scalar lowerRefineLevel_;
        scalar upperRefineLevel_;
        scalar unrefineLevel_;
        label maxRefinement_;
        label maxCells_;

        //- Number of cell layers around the interface band that are also
        //  refined (the refinement "safety margin"; distance ~ nBufferLayers x
        //  base-cell size). Read from dynamicT8codeFvMeshCoeffs::nBufferLayers.
        label nBufferLayers_;

        List<Pair<word>> correctFluxes_;

        //- Per base-cell refine flag (interface band dilated by nBufferLayers_);
        //  lifetime spans the adapt call that reads it
        std::vector<double> refineFlag_;

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

        //- Geometry of the collapsed (empty) direction for the 2D extrusion
        bool planar2D_;
        direction emptyDir_;
        scalar zMin_;
        scalar zMax_;

        //- Number of base (level-0) trees == number of base cells
        label nBaseTrees_;

        //- t8code objects (opaque; real types only in the .C)
        void* forest_;
        void* cmesh_;
        void* scheme_;
        bool  forestInitialized_;

        //- Per-element field storage (lifetime spans every adapt call)
        std::vector<double> alpha1Data_;
        std::vector<double> uData_;
        std::vector<double> pdData_;

        //- Field view passed to the C callbacks (points into the vectors above)
        T8FieldData fieldData_;


    // Private Member Functions

        void readDict();

        //- OF -> t8 : build the quad cmesh + level-0 forest from this mesh
        void initForest();

        //- OF -> t8 : copy alpha1/U/pd into the per-element vectors + user_data
        void copyFieldsToForest();

        //- Drive one t8code adapt (decision) + balance + field interpolation
        bool runT8codeAdapt();

        //- t8 -> OF : rebuild the WHOLE polyMesh from the forest leaves
        bool forestToOpenFOAM();

        //- t8 -> OF : rebuild only the cells t8 changed (+ finer-neighbour halo)
        bool forestToOpenFOAMIncremental();

        //- re-derive leafToCell_ (forest leaf <-> OF cell) by centroid hash
        void rebuildLeafToCell();

        //- t8 -> OF : write the interpolated per-element fields back
        void copyFieldsFromForest();

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

        explicit dynamicT8codeFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicT8codeFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


} // End namespace Foam

#endif

// ************************************************************************* //
