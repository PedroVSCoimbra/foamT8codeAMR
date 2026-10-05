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

\*---------------------------------------------------------------------------*/

#include "dynamicT8codeFvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceInterpolate.H"
#include "directTopoChange.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

#include <map>
#include <unordered_map>
#include <array>
#include <chrono>
#include <algorithm>

// t8code headers (kept out of the .H to avoid leaking C++17/MPI types)
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_helpers.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_vtk/t8_vtk.h>
#include <t8_schemes/t8_scheme.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>

#include <vector>

// Convenience casts from the opaque void* members to real t8code types
#define FOREST (static_cast<t8_forest_t>(forest_))
#define CMESH  (static_cast<t8_cmesh_t>(cmesh_))
#define SCHEME (static_cast<const t8_scheme_c*>(scheme_))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicT8codeFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicT8codeFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * //
//  t8code adapt callback (free C function)
//
//  user_data points to an AdaptCtx whose lifetime is owned by the caller
//  (runT8codeAdapt), guaranteeing it outlives this call.
// * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * //

namespace
{

struct AdaptCtx
{
    const double* alpha1;          // per OLD-element alpha1
    const double* refineFlag;      // per base-cell: >0.5 => refine (band+buffer)
    const Foam::label* treeOffset; // global leaf offset of each tree (old forest)
    double lower;
    double upper;
    int maxLevel;
};

int t8AdaptCallback
(
    t8_forest_t forest,
    t8_forest_t /*forestFrom*/,
    t8_locidx_t whichTree,
    const t8_eclass_t treeClass,
    t8_locidx_t lelementId,
    const t8_scheme_c* scheme,
    const int isFamily,
    const int numElements,
    t8_element_t* elements[]
)
{
    const AdaptCtx* ctx =
        static_cast<const AdaptCtx*>(t8_forest_get_user_data(forest));
    if (!ctx || !ctx->alpha1) return 0;

    const int level = t8_element_get_level(scheme, treeClass, elements[0]);
    const Foam::label gid = ctx->treeOffset[whichTree] + lelementId;

    // Per-leaf refine flag (refineFlag is indexed by global leaf == OF cell;
    // it already includes the interface band + nBufferLayers dilation). One
    // level of refinement per adapt call (recursive = 0); deeper levels build
    // up over successive refineIntervals, as in dynamicRefineFvMesh.
    const bool flagged = ctx->refineFlag && ctx->refineFlag[gid] > 0.5;

    // Coarsen a full family only if ALL siblings are clearly on the SAME side
    // of the interface (all gas, or all liquid). A family that straddles the
    // interface must NOT be coarsened: doing so destroys the interface and the
    // (injected) coarse value is wildly wrong -> blows up interFoam's
    // interface compression.
    if (isFamily && numElements > 1)
    {
        bool allLow = true;   // every sibling below the lower threshold (gas)
        bool allHigh = true;  // every sibling above the upper threshold (liquid)
        bool anyFlagged = false;
        for (int i = 0; i < numElements; ++i)
        {
            const double a = ctx->alpha1[gid + i];
            if (a >= ctx->lower) allLow = false;
            if (a <= ctx->upper) allHigh = false;
            if (ctx->refineFlag && ctx->refineFlag[gid + i] > 0.5)
            {
                anyFlagged = true;
            }
        }
        // Keep every family in the interface band and its buffer at its
        // current level. Coarsening it creates mixed levels along the
        // rectangular interface, especially where it reaches a wall.
        // Coarsen gas-only families; the liquid wake stays refined.
        // Coarsening pure-liquid families leaves impact zones under-
        // resolved and ejects numerical spray (dam-break wall run-up);
        // the native dynamicRefineFvMesh likewise preserves the wake.
        if (!anyFlagged && allLow && level > 0)
        {
            return -1;
        }
        // Refine when ANY member is flagged: t8 decides per family, so
        // testing only the first member leaves part of the band coarse
        // for a full adapt cycle, diffusing the interface (the missed
        // families sit inside the band at low level while neighbours
        // are already fine).
        if (anyFlagged && level < ctx->maxLevel)
        {
            return 1;
        }
        return 0;
    }

    // Refine (recursively) every element of a flagged cell up to maxLevel
    if (flagged && level < ctx->maxLevel)
    {
        return 1;
    }

    return 0;
}

// Cumulative leaf offset of every local tree -> size (numTrees+1)
Foam::labelList computeTreeOffsets(t8_forest_t forest)
{
    const t8_locidx_t nTrees = t8_forest_get_num_local_trees(forest);
    Foam::labelList off(nTrees + 1, 0);
    for (t8_locidx_t t = 0; t < nTrees; ++t)
    {
        off[t + 1] =
            off[t] + t8_forest_get_tree_num_leaf_elements(forest, t);
    }
    return off;
}

// Build the NEW-cell -> OLD-cell map straight from the t8 parent<->child
// relationship (no geometric search). OF cells are in forest-leaf order, so the
// old/new global leaf indices ARE the old/new OF cell indices.
//   refine 1->N : every child's master is the parent old cell
//   same   1->1 : master is the old cell
//   coarsen N->1: master is the first outgoing child (same-side, so equivalent)
struct ReplaceCtx
{
    Foam::label* master;            // size = new nLeaf
    const Foam::label* oldOff;
    const Foam::label* newOff;
    unsigned char* newDirty;        // size = new nLeaf (nullable)
    unsigned char* oldDirty;        // size = old nLeaf (nullable)

    const double* oldAlpha;
    const double* oldU;
    const double* oldPd;
    double* newAlpha;
    double* newU;
    double* newPd;
};

void t8ReplaceCallback
(
    t8_forest_t /*forest_old*/,
    t8_forest_t forest_new,
    t8_locidx_t whichTree,
    const t8_eclass_t /*tree_class*/,
    const t8_scheme_c* /*scheme*/,
    const int /*refine*/,
    const int numOutgoing,
    const t8_locidx_t firstOutgoing,
    const int numIncoming,
    const t8_locidx_t firstIncoming
)
{
    ReplaceCtx* c =
        static_cast<ReplaceCtx*>(t8_forest_get_user_data(forest_new));
    if (!c) return;

    const Foam::label go = c->oldOff[whichTree] + firstOutgoing;
    const Foam::label gn = c->newOff[whichTree] + firstIncoming;

    if (numIncoming >= numOutgoing)
    {
        for (int i = 0; i < numIncoming; ++i)
        {
            c->master[gn + i] = go;
            if (c->oldAlpha && c->newAlpha)
                c->newAlpha[gn + i] = c->oldAlpha[go];
            if (c->oldPd && c->newPd)
                c->newPd[gn + i] = c->oldPd[go];
            if (c->oldU && c->newU)
            {
                for (int cmpt = 0; cmpt < 3; ++cmpt)
                    c->newU[3*(gn + i) + cmpt] = c->oldU[3*go + cmpt];
            }
        }
    }
    else
    {
        c->master[gn] = go;
        if (c->oldAlpha && c->newAlpha)
        {
            double sum = 0.0;
            for (int i = 0; i < numOutgoing; ++i) sum += c->oldAlpha[go + i];
            c->newAlpha[gn] = sum/double(numOutgoing);
        }
        if (c->oldPd && c->newPd)
        {
            double sum = 0.0;
            for (int i = 0; i < numOutgoing; ++i) sum += c->oldPd[go + i];
            c->newPd[gn] = sum/double(numOutgoing);
        }
        if (c->oldU && c->newU)
        {
            for (int cmpt = 0; cmpt < 3; ++cmpt)
            {
                double sum = 0.0;
                for (int i = 0; i < numOutgoing; ++i)
                    sum += c->oldU[3*(go + i) + cmpt];
                c->newU[3*gn + cmpt] = sum/double(numOutgoing);
            }
        }
    }

    // dirty = refine (1->N) or coarsen (N->1); unchanged (1->1) stays clean
    const bool dirty = (numOutgoing != 1) || (numIncoming != 1);
    if (dirty)
    {
        if (c->newDirty)
            for (int i = 0; i < numIncoming; ++i) c->newDirty[gn + i] = 1;
        if (c->oldDirty)
            for (int i = 0; i < numOutgoing; ++i) c->oldDirty[go + i] = 1;
    }
}

// Hash for a quantised 3D point key (for the unordered_map point dedup)
struct Key3Hash
{
    std::size_t operator()(const std::array<long,3>& k) const noexcept
    {
        std::size_t h = 1469598103934665603ULL;
        for (int i = 0; i < 3; ++i)
        {
            h ^= std::size_t(k[i]);
            h *= 1099511628211ULL;
        }
        return h;
    }
};

} // End anonymous namespace


using namespace Foam;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicT8codeFvMesh::dynamicT8codeFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    refineInterval_(1),
    fieldName_("alpha1"),
    lowerRefineLevel_(0.001),
    upperRefineLevel_(0.999),
    unrefineLevel_(10),
    maxRefinement_(2),
    maxCells_(1000000),
    nBufferLayers_(1),
    correctFluxes_(),
    refineFlag_(),
    masterCell_(),
    incrementalRebuild_(false),
    leafToCell_(),
    leafDirty_(),
    oldDirty_(),
    curTimeIndex_(-1),
    planar2D_(false),
    emptyDir_(2),
    zMin_(0),
    zMax_(0),
    nBaseTrees_(0),
    forest_(nullptr),
    cmesh_(nullptr),
    scheme_(nullptr),
    forestInitialized_(false),
    alpha1Data_(),
    uData_(),
    pdData_(),
    fieldData_()
{
    // foam-extend initialises MPI via Pstream before the mesh is built; only
    // bootstrap it ourselves if (somehow) it is not yet up.
    int mpiInitialised = 0;
    MPI_Initialized(&mpiInitialised);
    if (!mpiInitialised)
    {
        int provided;
        sc_MPI_Init_thread
        (
            nullptr, nullptr, sc_MPI_THREAD_SINGLE, &provided
        );
    }

    readDict();

    Info<< "dynamicT8codeFvMesh: constructed on " << nCells()
        << " base cells" << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicT8codeFvMesh::~dynamicT8codeFvMesh()
{
    if (forest_)
    {
        // Clear any dangling user data before releasing the forest
        t8_forest_set_user_data(FOREST, nullptr);
        t8_forest_t f = FOREST;
        t8_forest_unref(&f);
        forest_ = nullptr;
    }
    if (scheme_)
    {
        const t8_scheme_c* s = SCHEME;
        t8_scheme_unref(const_cast<t8_scheme_c**>(&s));
        scheme_ = nullptr;
    }
    // cmesh_ is consumed by the forest (set_cmesh takes ownership); no destroy
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynamicT8codeFvMesh::readDict()
{
    IOdictionary dict
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");

    refineInterval_   = readLabel(coeffs.lookup("refineInterval"));
    fieldName_        = word(coeffs.lookup("field"));
    lowerRefineLevel_ = readScalar(coeffs.lookup("lowerRefineLevel"));
    upperRefineLevel_ = readScalar(coeffs.lookup("upperRefineLevel"));
    unrefineLevel_    = readScalar(coeffs.lookup("unrefineLevel"));
    maxRefinement_    = readLabel(coeffs.lookup("maxRefinement"));
    maxCells_         = readLabel(coeffs.lookup("maxCells"));
    nBufferLayers_    = coeffs.lookupOrDefault<label>("nBufferLayers", 1);
    incrementalRebuild_ =
        coeffs.lookupOrDefault<bool>("incrementalRebuild", false);
    correctFluxes_    = List<Pair<word>>(coeffs.lookup("correctFluxes"));
}


void dynamicT8codeFvMesh::initForest()
{
    nBaseTrees_ = nCells();

    planar2D_ = false;
    const polyBoundaryMesh& pbm = boundaryMesh();
    forAll(pbm, pI)
    {
        if (pbm[pI].type() == "empty")
        {
            planar2D_ = true;
            break;
        }
    }

    const boundBox bb(points());
    zMin_ = bb.min().z();
    zMax_ = bb.max().z();

    // ---- Build the cmesh: QUAD for empty 2D, HEX for full 3D --------------
    t8_cmesh_t cmesh;
    t8_cmesh_init(&cmesh);
    t8_cmesh_register_geometry(cmesh, t8_geometry_linear_new());

    const pointField& pts = points();

    forAll(cells(), cellI)
    {
        const labelList cPoints = cells()[cellI].labels(faces());

        // cell centroid, used to place each corner in t8 hex z-order
        point cc(0, 0, 0);
        forAll(cPoints, k) cc += pts[cPoints[k]];
        cc /= max(cPoints.size(), 1);

        if (planar2D_)
        {
            // t8 quad corner index = (x > cx) + 2*(y > cy).  Keep z fixed in
            // t8code and extrude the resulting quad leaves back to one hex
            // layer during OpenFOAM reconstruction.
            point v[4];
            forAll(cPoints, k)
            {
                const point& p = pts[cPoints[k]];
                const int idx =
                    (p.x() > cc.x() ? 1 : 0)
                  + (p.y() > cc.y() ? 2 : 0);
                v[idx] = point(p.x(), p.y(), 0.0);
            }

            double treeVerts[12];
            for (int corner = 0; corner < 4; ++corner)
            {
                treeVerts[3*corner + 0] = v[corner].x();
                treeVerts[3*corner + 1] = v[corner].y();
                treeVerts[3*corner + 2] = v[corner].z();
            }

            t8_cmesh_set_tree_class(cmesh, cellI, T8_ECLASS_QUAD);
            t8_cmesh_set_tree_vertices(cmesh, cellI, treeVerts, 4);
        }
        else
        {
            // t8 hex corner index = (x>cx) + 2*(y>cy) + 4*(z>cz)
            point v[8];
            forAll(cPoints, k)
            {
                const point& p = pts[cPoints[k]];
                const int idx =
                    (p.x() > cc.x() ? 1 : 0)
                  + (p.y() > cc.y() ? 2 : 0)
                  + (p.z() > cc.z() ? 4 : 0);
                v[idx] = p;
            }

            double treeVerts[24];
            for (int corner = 0; corner < 8; ++corner)
            {
                treeVerts[3*corner + 0] = v[corner].x();
                treeVerts[3*corner + 1] = v[corner].y();
                treeVerts[3*corner + 2] = v[corner].z();
            }

            t8_cmesh_set_tree_class(cmesh, cellI, T8_ECLASS_HEX);
            t8_cmesh_set_tree_vertices(cmesh, cellI, treeVerts, 8);
        }
    }

    // Auto-compute the face joins from the tree vertices in the stash
    t8_cmesh_set_join_by_stash(cmesh, nullptr, 1);

    t8_cmesh_commit(cmesh, sc_MPI_COMM_WORLD);
    cmesh_ = cmesh;

    // ---- Build the level-0 forest -----------------------------------------
    scheme_ = const_cast<t8_scheme_c*>(t8_scheme_new_default());

    t8_forest_t forest;
    t8_forest_init(&forest);
    t8_forest_set_cmesh(forest, cmesh, sc_MPI_COMM_WORLD);
    t8_forest_set_scheme(forest, SCHEME);
    t8_forest_set_level(forest, 0);
    t8_forest_commit(forest);
    forest_ = forest;

    forestInitialized_ = true;

    Info<< "dynamicT8codeFvMesh::initForest: " << nBaseTrees_
        << (planar2D_ ? " quad trees (2D empty)" : " hex trees (3D)")
        << endl;

    // ---- read the mesh back THROUGH t8 and print what t8 sees -------------
    {
        const t8_scheme_c* sch = SCHEME;
        const t8_locidx_t nTrees =
            t8_forest_get_num_local_trees(FOREST);
        const t8_locidx_t nElem =
            t8_forest_get_local_num_leaf_elements(FOREST);

        Info<< nl
            << "================ T8CODE MESH VIEW =================" << nl
            << "  local trees (= base cells) : " << label(nTrees) << nl
            << "  leaf elements              : " << label(nElem) << nl
            << "  (OF reports cells          : " << nCells() << ")" << nl;

        // accumulate t8's idea of the geometry while walking the leaves
        scalar t8Vol = 0;
        point bbMin( GREAT,  GREAT,  GREAT);
        point bbMax(-GREAT, -GREAT, -GREAT);
        label shown = 0;

        for (t8_locidx_t tloc = 0; tloc < nTrees; ++tloc)
        {
            const t8_locidx_t nE =
                t8_forest_get_tree_num_leaf_elements(FOREST, tloc);
            for (t8_locidx_t e = 0; e < nE; ++e)
            {
                const t8_element_t* el =
                    t8_forest_get_leaf_element_in_tree(FOREST, tloc, e);
                    const t8_eclass_t treeClass =
                        (planar2D_ ? T8_ECLASS_QUAD : T8_ECLASS_HEX);
                    const int lvl =
                        t8_element_get_level(sch, treeClass, el);

                point c[8];
                const int nCorners = (planar2D_ ? 4 : 8);
                for (int corner = 0; corner < nCorners; ++corner)
                {
                    double xyz[3];
                    t8_forest_element_coordinate(FOREST, tloc, el, corner, xyz);
                    c[corner] = point(xyz[0], xyz[1], xyz[2]);
                    bbMin = min(bbMin, c[corner]);
                    bbMax = max(bbMax, c[corner]);
                }
                const scalar dx = mag(c[1] - c[0]);
                const scalar dy = mag(c[2] - c[0]);
                const scalar dz = planar2D_ ? (zMax_ - zMin_) : mag(c[4] - c[0]);
                t8Vol += dx*dy*dz;

                if (shown < 4)
                {
                    Info<< "  element " << label(shown)
                        << " : tree " << label(tloc)
                        << ", level " << lvl
                        << ", v0 " << c[0] << ", dx " << dx
                        << ", dy " << dy << ", dz " << dz << nl;
                    ++shown;
                }
            }
        }

        Info<< "  bounding box (from t8)     : " << bbMin
            << " -> " << bbMax << nl
            << "  t8 total volume            : " << t8Vol
            << "  (OF total volume = " << gSum(V()) << ")" << nl
            << "==================================================" << nl
            << endl;
    }
}


void dynamicT8codeFvMesh::copyFieldsToForest()
{
    const t8_locidx_t nElem = t8_forest_get_local_num_leaf_elements(FOREST);

    alpha1Data_.assign(nElem, 0.0);
    uData_.assign(std::size_t(nElem) * 3, 0.0);
    pdData_.assign(nElem, 0.0);

    // Map forest leaf -> OF cell through leafToCell_. On the very first adapt
    // (level-0 forest, no rebuild yet) leafToCell_ is empty -> identity.
    const bool useMap = (label(leafToCell_.size()) == nElem);
    auto cellOf = [&](label leaf) -> label
    { return useMap ? leafToCell_[leaf] : leaf; };

    const volScalarField& a = lookupObject<volScalarField>(fieldName_);
    for (label i = 0; i < nElem; ++i)
    {
        const label c = cellOf(i);
        if (c >= 0 && c < a.size()) alpha1Data_[i] = a[c];
    }

    if (foundObject<volVectorField>("U"))
    {
        const volVectorField& U = lookupObject<volVectorField>("U");
        for (label i = 0; i < nElem; ++i)
        {
            const label c = cellOf(i);
            if (c >= 0 && c < U.size())
            {
                uData_[3*i + 0] = U[c].x();
                uData_[3*i + 1] = U[c].y();
                uData_[3*i + 2] = U[c].z();
            }
        }
    }
    if (foundObject<volScalarField>("pd"))
    {
        const volScalarField& pd = lookupObject<volScalarField>("pd");
        for (label i = 0; i < nElem; ++i)
        {
            const label c = cellOf(i);
            if (c >= 0 && c < pd.size()) pdData_[i] = pd[c];
        }
    }

    fieldData_.alpha1 = alpha1Data_.data();
    fieldData_.U      = uData_.data();
    fieldData_.pd     = pdData_.data();
    fieldData_.nElements = int(nElem);
}


bool dynamicT8codeFvMesh::runT8codeAdapt()
{
    copyFieldsToForest();

    t8_forest_t forestOld = FOREST;

    const labelList oldTreeOffset = computeTreeOffsets(forestOld);
    const t8_locidx_t nOld =
        t8_forest_get_local_num_leaf_elements(forestOld);

    // ---- refine flag: interface band dilated by nBufferLayers ---------------
    // (the "safety margin": distance ~ nBufferLayers x cell size). Computed in
    // OF-CELL space (band on the OF field, dilation over OF faces), then mapped
    // to per-LEAF (what the adapt callback indexes) via leafToCell_. After an
    // incremental rebuild OF cell index != forest leaf index, so the map is
    // essential; on a from-scratch mesh leafToCell_ is the identity.
    {
        const label nc = nCells();
        std::vector<double> cellFlag(nc, 0.0);
        const volScalarField& a = lookupObject<volScalarField>(fieldName_);

        // Match dynamicRefineFvMesh exactly: interpolate cell values to
        // points, calculate the point error, and take the maximum error over
        // each cell's points.
        scalarField pointAlpha(points().size(), 0.0);
        labelList pointCount(points().size(), 0);
        const labelListList& pCells = pointCells();
        forAll(pCells, pointI)
        {
            forAll(pCells[pointI], i)
            {
                pointAlpha[pointI] += a[pCells[pointI][i]];
                pointCount[pointI]++;
            }
            if (pointCount[pointI] > 0)
                pointAlpha[pointI] /= pointCount[pointI];
        }

        const scalar halfLevel =
            0.5*(lowerRefineLevel_ + upperRefineLevel_);
        scalarField pointError(points().size(), -1.0);
        forAll(pointError, pointI)
        {
            if
            (
                pointAlpha[pointI] >= lowerRefineLevel_
             && pointAlpha[pointI] < upperRefineLevel_
            )
            {
                pointError[pointI] = mag(pointAlpha[pointI] - halfLevel);
            }
        }

        const labelListList& cPoints = cellPoints();
        forAll(cPoints, c)
        {
            scalar cellError = -1.0;
            forAll(cPoints[c], i)
            {
                cellError = max(cellError, pointError[cPoints[c][i]]);
            }
            if (cellError > 0.0)
                cellFlag[c] = 1.0;
        }

        // Boundary faces have no neighbour cell for the usual face-based
        // buffer. Protect their owner cells explicitly when an interface point
        // lies on the patch, while allowing protection to disappear when the
        // interface moves away on a later adaptation.
        const labelUList& fOwn = faceOwner();
        for (label faceI = nInternalFaces(); faceI < nFaces(); ++faceI)
        {
            bool interfaceOnPatch = false;
            const face& boundaryFace = faces()[faceI];
            forAll(boundaryFace, i)
            {
                if (pointError[boundaryFace[i]] > 0.0)
                {
                    interfaceOnPatch = true;
                    break;
                }
            }
            if (interfaceOnPatch)
                cellFlag[fOwn[faceI]] = 1.0;
        }

        const labelUList& fNei = faceNeighbour();
        for (label L = 0; L < nBufferLayers_; ++L)
        {
            std::vector<double> nf = cellFlag;
            forAll(fNei, f) // internal faces only
            {
                if (cellFlag[fOwn[f]] > 0.5) nf[fNei[f]] = 1.0;
                if (cellFlag[fNei[f]] > 0.5) nf[fOwn[f]] = 1.0;
            }
            cellFlag.swap(nf);
        }

        refineFlag_.assign(nOld, 0.0);
        const bool useMap = (label(leafToCell_.size()) == nOld);
        for (label leaf = 0; leaf < nOld; ++leaf)
        {
            const label c = useMap ? leafToCell_[leaf] : leaf;
            if (c >= 0 && c < nc) refineFlag_[leaf] = cellFlag[c];
        }
    }

    AdaptCtx ctx;
    ctx.alpha1     = alpha1Data_.data();
    ctx.refineFlag = refineFlag_.data();
    ctx.treeOffset = oldTreeOffset.begin();
    ctx.lower      = lowerRefineLevel_;
    ctx.upper      = upperRefineLevel_;
    ctx.maxLevel   = int(maxRefinement_);

    // Keep the old forest alive: t8_forest_new_adapt consumes one reference of
    // its argument. We also need the old forest for the parent<->child map.
    t8_forest_ref(forestOld);
    t8_forest_set_user_data(forestOld, &ctx);

    // 1) Adapt one level per call, as dynamicRefineFvMesh does. The next
    //    update recomputes the interface on the newly mapped mesh.
    t8_forest_t forestAdapt =
        t8_forest_new_adapt(forestOld, t8AdaptCallback, 0, 0, &ctx);

    // 2) 2:1 balance + face-ghost layer (needed by leaf_face_neighbors).
    //    NOTE: the balance is the dominant cost of the whole adapt (~400 ms at
    //    ~80k cells); it is mandatory (skipping it leaves >2:1 jumps that the
    //    reconstruction cannot represent -> open cells). t8 offers no public
    //    incremental balance, so this is a fixed floor comparable to hexRef8's
    //    entire mesh.update.
    t8_forest_t forestNew;
    t8_forest_init(&forestNew);
    t8_forest_set_balance(forestNew, forestAdapt, 0); // consumes forestAdapt
    t8_forest_set_ghost(forestNew, 1, T8_GHOST_FACES);
    t8_forest_commit(forestNew);

    const t8_locidx_t nNew =
        t8_forest_get_local_num_leaf_elements(forestNew);

    // 3) Build new-cell -> old-cell map from the t8 parent<->child relation
    //    (no geometric findCell). Used as the directTopoChange masterCell.
    const labelList newTreeOffset = computeTreeOffsets(forestNew);
    masterCell_.assign(nNew, -1);
    leafDirty_.assign(nNew, 0);
    oldDirty_.assign(nOld, 0);
    std::vector<double> newAlpha(nNew, 0.0);
    std::vector<double> newU(std::size_t(nNew)*3, 0.0);
    std::vector<double> newPd(nNew, 0.0);
    ReplaceCtx rctx;
    rctx.master = masterCell_.data();
    rctx.oldOff = oldTreeOffset.begin();
    rctx.newOff = newTreeOffset.begin();
    rctx.newDirty = leafDirty_.data();
    rctx.oldDirty = oldDirty_.data();
    rctx.oldAlpha = alpha1Data_.data();
    rctx.oldU = uData_.data();
    rctx.oldPd = pdData_.data();
    rctx.newAlpha = newAlpha.data();
    rctx.newU = newU.data();
    rctx.newPd = newPd.data();
    t8_forest_set_user_data(forestNew, &rctx);
    t8_forest_iterate_replace(forestNew, forestOld, t8ReplaceCallback);
    t8_forest_set_user_data(forestNew, nullptr);

    alpha1Data_.swap(newAlpha);
    uData_.swap(newU);
    pdData_.swap(newPd);

    t8_forest_set_user_data(forestOld, nullptr);
    t8_forest_unref(&forestOld);   // release the extra ref taken above
    forest_ = forestNew;

    Info<< "t8code adapt: " << label(nOld) << " -> " << label(nNew)
        << " elements" << endl;

    // Report a change whenever the layout changed, not only when the
    // element count changed. Same-count adapts (refine here, coarsen
    // there) still advance forest_ : skipping the OF rebuild then would
    // leave the mesh and leafToCell_ behind the forest, silently
    // scrambling every later field transfer (spray + mass error).
    bool layoutChanged = (nNew != nOld);
    for (label i = 0; i < nNew && !layoutChanged; ++i)
    {
        if (leafDirty_[std::size_t(i)]) layoutChanged = true;
    }
    for (label i = 0; i < nOld && !layoutChanged; ++i)
    {
        if (oldDirty_[std::size_t(i)]) layoutChanged = true;
    }
    return layoutChanged;
}


void dynamicT8codeFvMesh::correctMappedFluxes(const mapPolyMesh& map)
{
    // Mirror dynamicRefineFvMesh::update: the mapped flux (limited by MULES
    // before the topology change) is kept wherever an old master exists;
    // only genuinely new faces are recomputed from the mapped velocity.
    const labelList& faceMap = map.faceMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    forAll(correctFluxes_, i)
    {
        const word& phiName = correctFluxes_[i][0];
        const word& Uname   = correctFluxes_[i][1];
        if
        (
            foundObject<surfaceScalarField>(phiName)
         && foundObject<volVectorField>(Uname)
        )
        {
            surfaceScalarField& phi = const_cast<surfaceScalarField&>
                (lookupObject<surfaceScalarField>(phiName));
            const volVectorField& U =
                lookupObject<volVectorField>(Uname);
            const surfaceScalarField phiU((fvc::interpolate(U) & Sf()));

            auto needsRecompute = [&](label faceI) -> bool
            {
                const label oldFaceI =
                    faceI < faceMap.size() ? faceMap[faceI] : -1;
                return
                (
                    oldFaceI < 0
                 || oldFaceI >= reverseFaceMap.size()
                 || reverseFaceMap[oldFaceI] != faceI
                );
            };

            // Internal faces: recompute only unmapped/extra-child faces
            for (label faceI = 0; faceI < nInternalFaces(); ++faceI)
            {
                if (needsRecompute(faceI))
                {
                    phi[faceI] = phiU[faceI];
                }
            }

            // Boundary faces: same rule
            forAll(phi.boundaryField(), patchI)
            {
                fvsPatchScalarField& patchPhi =
                    phi.boundaryField()[patchI];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchI];
                label faceI = patchPhi.patch().patch().start();
                forAll(patchPhi, k)
                {
                    if (needsRecompute(faceI))
                    {
                        patchPhi[k] = patchPhiU[k];
                    }
                    ++faceI;
                }
            }
        }
    }
}


bool dynamicT8codeFvMesh::forestToOpenFOAM()
{
    t8_forest_t forest = FOREST;
    const t8_scheme_c* scheme = SCHEME;

    const t8_locidx_t nTrees = t8_forest_get_num_local_trees(forest);
    const t8_locidx_t nLeaf  = t8_forest_get_local_num_leaf_elements(forest);
    const labelList off = computeTreeOffsets(forest);

    if (planar2D_)
    {
        const t8_eclass_t treeClass = T8_ECLASS_QUAD;

        const boundBox& bb = bounds();
        const point bMin = bb.min();
        const point bMax = bb.max();

        auto sideOf = [&](const point& p) -> int
        {
            const scalar d[6] =
            {
                mag(p.x()-bMin.x()), mag(p.x()-bMax.x()),
                mag(p.y()-bMin.y()), mag(p.y()-bMax.y()),
                mag(p.z()-zMin_), mag(p.z()-zMax_)
            };
            int s = 0;
            for (int i = 1; i < 6; ++i) if (d[i] < d[s]) s = i;
            return s;
        };

        const polyBoundaryMesh& pbm = boundaryMesh();
        // Map each domain side to its patch by inspecting every boundary
        // face centre. Averaging a patch (e.g. empty covering front+back)
        // collapses opposite sides to the domain centre and leaves one side
        // unmapped, sending its faces to patch 0.
        labelList sidePatchID(6, -1);
        forAll(pbm, pI)
        {
            const vectorField& fc = pbm[pI].faceCentres();
            forAll(fc, k) sidePatchID[sideOf(fc[k])] = pI;
        }
        auto assignSide = [&](const point& fc) -> label
        {
            const label p = sidePatchID[sideOf(fc)];
            return p >= 0 ? p : 0;
        };

        directTopoChange meshMod(pbm.size());

        // OLD-mesh face data (valid until changeMesh below). Each new
        // face maps from the old face with the closest centre among its
        // master cells' faces, so mapped surface fields stay finite.
        // Without a master every faceMap entry is -1 and mapped surface
        // fields (e.g. nHatf_) are left uninitialized, producing NaN and
        // a SIGFPE in the flux weights.
        const cellList& oldCellFaces = cells();
        vectorField oldFC(nFaces());
        {
            const vectorField& ifc = Cf().internalField();
            forAll(ifc, f) oldFC[f] = ifc[f];
            const polyBoundaryMesh& obm = boundaryMesh();
            forAll(obm, p)
            {
                const vectorField& pfc = obm[p].faceCentres();
                const label s = obm[p].start();
                forAll(pfc, i) oldFC[s + i] = pfc[i];
            }
        }

        std::unordered_map<std::array<long,3>, label, Key3Hash> pointMap;
        pointMap.reserve(std::size_t(nLeaf) * 6);
        const double qtol = 1e-9;
        auto getPoint = [&](const point& p) -> label
        {
            std::array<long,3> key =
            {
                std::lround(p.x()/qtol),
                std::lround(p.y()/qtol),
                std::lround(p.z()/qtol)
            };
            auto it = pointMap.find(key);
            if (it != pointMap.end()) return it->second;
            const label id = meshMod.addPoint(p, -1, -1, true);
            pointMap[key] = id;
            return id;
        };

        auto addPolyFace = [&]
        (
            std::vector<point> poly,
            const label own, const label nei, const label patch,
            const point& ownC, const point& orientRef, const bool boundary
        )
        {
            const int np = int(poly.size());
            if (np < 3) return;

            point mn = poly[0], mx = poly[0];
            for (const point& p : poly) { mn = min(mn, p); mx = max(mx, p); }
            const vector ext = mx - mn;
            int axis = 0;
            if (ext[1] < ext[axis]) axis = 1;
            if (ext[2] < ext[axis]) axis = 2;
            const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;

            point ctr(0,0,0);
            for (const point& p : poly) ctr += p;
            ctr /= scalar(np);

            std::sort(poly.begin(), poly.end(), [&](const point& A, const point& B)
            {
                return std::atan2(A[a2]-ctr[a2], A[a1]-ctr[a1])
                     < std::atan2(B[a2]-ctr[a2], B[a1]-ctr[a1]);
            });

            vector nrm(0,0,0);
            for (int i = 0; i < np; ++i)
            {
                const point& P = poly[i];
                const point& Q = poly[(i+1)%np];
                nrm.x() += (P.y()-Q.y())*(P.z()+Q.z());
                nrm.y() += (P.z()-Q.z())*(P.x()+Q.x());
                nrm.z() += (P.x()-Q.x())*(P.y()+Q.y());
            }
            const vector dir = boundary ? (ctr - ownC) : (orientRef - ownC);
            const bool rev = ((nrm & dir) < 0);

            // Face master: old face with the closest centre among the
            // master cells' faces. New cells are forest leaves, so the
            // old cell comes from masterCell_ via leafToCell_.
            label masterFace = -1;
            {
                const label nwc[2] = {own, nei};
                label mc[2] = {-1, -1};
                for (int k = 0; k < 2; ++k)
                {
                    const label nc = nwc[k];
                    if (nc >= 0 && nc < label(masterCell_.size()))
                    {
                        const label ml = masterCell_[nc];
                        // Same convention as the cell masters below:
                        // fall back to the leaf index itself when the
                        // leaf->cell map misses (first rebuild).
                        label oc = -1;
                        if
                        (
                            ml >= 0
                         && ml < label(leafToCell_.size())
                        )
                        {
                            oc = leafToCell_[ml];
                        }
                        else
                        {
                            oc = ml;
                        }
                        if (oc >= 0 && oc < oldCellFaces.size())
                        {
                            mc[k] = oc;
                        }
                    }
                }
                scalar bestD = GREAT;
                for (int k = 0; k < 2; ++k)
                {
                    if (mc[k] < 0) continue;
                    if (k == 1 && mc[1] == mc[0]) continue;
                    const cell& cFaces = oldCellFaces[mc[k]];
                    forAll(cFaces, fi)
                    {
                        const label of = cFaces[fi];
                        if (of < 0 || of >= oldFC.size()) continue;
                        const scalar d = mag(oldFC[of] - ctr);
                        if (d < bestD)
                        {
                            bestD = d;
                            masterFace = of;
                        }
                    }
                }
            }

            labelList fp(np);
            for (int i = 0; i < np; ++i)
                fp[i] = getPoint(rev ? poly[np-1-i] : poly[i]);
            meshMod.addFace
            (
                face(fp), own, nei, -1, -1, masterFace,
                false, patch, -1, false
            );
        };

        std::vector<point> leafCentre(nLeaf);
        std::vector<std::array<point,8>> leafCorners(nLeaf);
        for (t8_locidx_t t = 0; t < nTrees; ++t)
        {
            const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
            for (t8_locidx_t e = 0; e < nE; ++e)
            {
                const label gc = off[t] + e;
                const t8_element_t* leaf =
                    t8_forest_get_leaf_element_in_tree(forest, t, e);
                point c(0,0,0);
                for (int corner = 0; corner < 4; ++corner)
                {
                    double xyz[3];
                    t8_forest_element_coordinate(forest, t, leaf, corner, xyz);
                    const point p0(xyz[0], xyz[1], zMin_);
                    const point p1(xyz[0], xyz[1], zMax_);
                    leafCorners[gc][corner] = p0;
                    leafCorners[gc][corner + 4] = p1;
                    c += p0 + p1;
                    getPoint(p0);
                    getPoint(p1);
                }
                c /= 8.0;
                leafCentre[gc] = c;

                label master = -1;
                if (gc < label(masterCell_.size()))
                {
                    const label masterLeaf = masterCell_[gc];
                    if
                    (
                        masterLeaf >= 0
                     && masterLeaf < label(leafToCell_.size())
                    )
                    {
                        master = leafToCell_[masterLeaf];
                    }
                    else
                    {
                        master = masterLeaf;
                    }
                }
                meshMod.addCell(-1, -1, -1, master, -1);
            }
        }

        for (t8_locidx_t t = 0; t < nTrees; ++t)
        {
            const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
            for (t8_locidx_t e = 0; e < nE; ++e)
            {
                const label gc = off[t] + e;
                const t8_element_t* leaf =
                    t8_forest_get_leaf_element_in_tree(forest, t, e);
                const int myLevel = t8_element_get_level(scheme, treeClass, leaf);
                const std::array<point,8>& cc = leafCorners[gc];

                for (int k = 0; k < 4; ++k)
                {
                    int c0 = t8_element_get_face_corner(scheme, treeClass, leaf, k, 0);
                    int c1 = t8_element_get_face_corner(scheme, treeClass, leaf, k, 1);

                    t8_element_t** neighs = nullptr;
                    int* df = nullptr; int nN = 0;
                    t8_locidx_t* ei = nullptr; t8_eclass_t ec;
                    t8_forest_leaf_face_neighbors
                    (
                        forest, t, leaf,
                        const_cast<const t8_element_t***>(&neighs),
                        k, &df, &nN, &ei, &ec
                    );

                    bool emit = false, boundary = false;
                    label own = gc, nei = -1;
                    if (nN == 0)
                    {
                        emit = true;
                        boundary = true;
                    }
                    else if (nN == 1)
                    {
                        const label ng = ei[0];
                        const int nl = t8_element_get_level(scheme, ec, neighs[0]);
                        if (ng < nLeaf && (nl < myLevel || (nl == myLevel && gc < ng)))
                        {
                            emit = true;
                            own = min(gc, ng);
                            nei = max(gc, ng);
                        }
                    }
                    if (nN > 0) { T8_FREE(neighs); T8_FREE(ei); T8_FREE(df); }
                    if (!emit) continue;

                    std::vector<point> poly;
                    poly.push_back(cc[c0]);
                    poly.push_back(cc[c1]);
                    poly.push_back(cc[c1 + 4]);
                    poly.push_back(cc[c0 + 4]);
                    point fcen(0,0,0);
                    for (const point& p : poly) fcen += p;
                    fcen /= scalar(poly.size());
                    const label patch = boundary ? assignSide(fcen) : -1;
                    addPolyFace(poly, own, nei, patch, leafCentre[own],
                        boundary ? leafCentre[own] : leafCentre[nei], boundary);
                }

                const int front[4] = {0, 1, 3, 2};
                const int back[4] = {4, 6, 7, 5};
                for (int side = 0; side < 2; ++side)
                {
                    std::vector<point> poly;
                    const int* ids = (side == 0 ? front : back);
                    for (int i = 0; i < 4; ++i) poly.push_back(cc[ids[i]]);
                    point fcen(0,0,0);
                    for (const point& p : poly) fcen += p;
                    fcen /= scalar(poly.size());
                    addPolyFace(poly, gc, -1, assignSide(fcen), leafCentre[gc],
                        leafCentre[gc], true);
                }
            }
        }

        autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);
        updateMesh(map);
        rebuildLeafToCell();

        {
            HashTable<volScalarField*> sFlds =
                objectRegistry::lookupClass<volScalarField>();
            forAllIter(HashTable<volScalarField*>, sFlds, it)
            {
                const word& name = it.key();
                if (name != "nut")
                {
                    (*it)->correctBoundaryConditions();
                }
            }
            HashTable<volVectorField*> vFlds =
                objectRegistry::lookupClass<volVectorField>();
            forAllIter(HashTable<volVectorField*>, vFlds, it)
                (*it)->correctBoundaryConditions();
        }

        // Keep the mapped (MULES-limited) fluxes; recompute from the
        // mapped velocity only faces without an old master (native rule).
        correctMappedFluxes(map);

        Info<< "dynamicT8codeFvMesh: rebuilt 2D empty OF mesh with "
            << nCells() << " cells, " << nFaces() << " faces, "
            << nPoints() << " points" << endl;

        return true;
    }

    // ---- map each of the 6 domain sides to its boundary patch -------------
    const boundBox& bb = bounds();
    const point bMin = bb.min();
    const point bMax = bb.max();

    auto sideOf = [&](const point& p) -> int
    {
        const scalar d[6] =
        {
            mag(p.x()-bMin.x()), mag(p.x()-bMax.x()),
            mag(p.y()-bMin.y()), mag(p.y()-bMax.y()),
            mag(p.z()-bMin.z()), mag(p.z()-bMax.z())
        };
        int s = 0;
        for (int i = 1; i < 6; ++i) if (d[i] < d[s]) s = i;
        return s;
    };

    const polyBoundaryMesh& pbm = boundaryMesh();
    // Per-face mapping (see 2D branch): averaging collapses patches that
    // span opposite sides (e.g. empty front+back) to the centre.
    labelList sidePatchID(6, -1);
    forAll(pbm, pI)
    {
        const vectorField& fc = pbm[pI].faceCentres();
        forAll(fc, k) sidePatchID[sideOf(fc[k])] = pI;
    }
    auto assignSide = [&](const point& fc) -> label
    {
        const label p = sidePatchID[sideOf(fc)];
        return p >= 0 ? p : 0;
    };

    // ---- topology builder (no starting mesh) ------------------------------
    directTopoChange meshMod(pbm.size());

    // OLD-mesh face data (valid until changeMesh below). Each new face
    // maps from the old face with the closest centre among its master
    // cells' faces, so mapped surface fields stay finite (see 2D
    // branch: unmapped faces leave e.g. nHatf_ uninitialized -> NaN).
    const cellList& oldCellFaces = cells();
    vectorField oldFC(nFaces());
    {
        const vectorField& ifc = Cf().internalField();
        forAll(ifc, f) oldFC[f] = ifc[f];
        const polyBoundaryMesh& obm = boundaryMesh();
        forAll(obm, p)
        {
            const vectorField& pfc = obm[p].faceCentres();
            const label s = obm[p].start();
            forAll(pfc, i) oldFC[s + i] = pfc[i];
        }
    }

    std::unordered_map<std::array<long,3>, label, Key3Hash> pointMap;
    pointMap.reserve(std::size_t(nLeaf) * 4);
    const double qtol = 1e-9;
    auto getPoint = [&](const point& p) -> label
    {
        std::array<long,3> key =
        {
            std::lround(p.x()/qtol),
            std::lround(p.y()/qtol),
            std::lround(p.z()/qtol)
        };
        auto it = pointMap.find(key);
        if (it != pointMap.end()) return it->second;
        const label id = meshMod.addPoint(p, -1, -1, true);
        pointMap[key] = id;
        return id;
    };

    // does a point already exist? (used to detect 3D hanging nodes: a coarse
    // edge midpoint exists iff a finer neighbour cell created it as a corner)
    auto pointExists = [&](const point& p) -> bool
    {
        const std::array<long,3> key =
        {
            std::lround(p.x()/qtol),
            std::lround(p.y()/qtol),
            std::lround(p.z()/qtol)
        };
        return pointMap.find(key) != pointMap.end();
    };

    // add a planar (axis-aligned) polygon face, ordered + oriented correctly
    auto addPolyFace = [&]
    (
        std::vector<point> poly,            // by value: we sort it
        const label own, const label nei, const label patch,
        const point& ownC, const point& orientRef, const bool boundary
    )
    {
        const int np = int(poly.size());
        if (np < 3) return;

        // normal axis = the coordinate with (near) zero extent over the face
        point mn = poly[0], mx = poly[0];
        for (const point& p : poly) { mn = min(mn, p); mx = max(mx, p); }
        const vector ext = mx - mn;
        int axis = 0;
        if (ext[1] < ext[axis]) axis = 1;
        if (ext[2] < ext[axis]) axis = 2;
        const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;

        point ctr(0,0,0);
        for (const point& p : poly) ctr += p;
        ctr /= scalar(np);

        std::sort
        (
            poly.begin(), poly.end(),
            [&](const point& A, const point& B)
            {
                return std::atan2(A[a2]-ctr[a2], A[a1]-ctr[a1])
                     < std::atan2(B[a2]-ctr[a2], B[a1]-ctr[a1]);
            }
        );

        // Newell normal of the now-ordered loop
        vector nrm(0,0,0);
        for (int i = 0; i < np; ++i)
        {
            const point& P = poly[i];
            const point& Q = poly[(i+1)%np];
            nrm.x() += (P.y()-Q.y())*(P.z()+Q.z());
            nrm.y() += (P.z()-Q.z())*(P.x()+Q.x());
            nrm.z() += (P.x()-Q.x())*(P.y()+Q.y());
        }
        const vector dir = boundary ? (ctr - ownC) : (orientRef - ownC);
        const bool rev = ((nrm & dir) < 0);

        // Face master: old face with the closest centre among the
        // master cells' faces (same convention as the 2D branch).
        label masterFace = -1;
        {
            const label nwc[2] = {own, nei};
            label mc[2] = {-1, -1};
            for (int k = 0; k < 2; ++k)
            {
                const label nc = nwc[k];
                if (nc >= 0 && nc < label(masterCell_.size()))
                {
                    const label ml = masterCell_[nc];
                    label oc = -1;
                    if
                    (
                        ml >= 0
                     && ml < label(leafToCell_.size())
                    )
                    {
                        oc = leafToCell_[ml];
                    }
                    else
                    {
                        oc = ml;
                    }
                    if (oc >= 0 && oc < oldCellFaces.size())
                    {
                        mc[k] = oc;
                    }
                }
            }
            scalar bestD = GREAT;
            for (int k = 0; k < 2; ++k)
            {
                if (mc[k] < 0) continue;
                if (k == 1 && mc[1] == mc[0]) continue;
                const cell& cFaces = oldCellFaces[mc[k]];
                forAll(cFaces, fi)
                {
                    const label of = cFaces[fi];
                    if (of < 0 || of >= oldFC.size()) continue;
                    const scalar d = mag(oldFC[of] - ctr);
                    if (d < bestD)
                    {
                        bestD = d;
                        masterFace = of;
                    }
                }
            }
        }

        labelList fp(np);
        for (int i = 0; i < np; ++i)
            fp[i] = getPoint(rev ? poly[np-1-i] : poly[i]);
        meshMod.addFace
        (
            face(fp), own, nei, -1, -1, masterFace,
            false, patch, -1, false
        );
    };

    // PROFILE: sub-phase timers
    auto _pt = std::chrono::steady_clock::now();
    auto _ms = [&]()
    {
        auto n = std::chrono::steady_clock::now();
        double d = std::chrono::duration<double,std::milli>(n - _pt).count();
        _pt = n; return d;
    };

    // ---- pass 1: cells + centroids ----------------------------------------
    std::vector<point> leafCentre(nLeaf);
    std::vector<std::array<point,8>> leafCorners(nLeaf); // cache (reused pass 2)
    for (t8_locidx_t t = 0; t < nTrees; ++t)
    {
        const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
        for (t8_locidx_t e = 0; e < nE; ++e)
        {
            const label gc = off[t] + e;
            const t8_element_t* leaf =
                t8_forest_get_leaf_element_in_tree(forest, t, e);
            point c(0,0,0);
            for (int corner = 0; corner < 8; ++corner)
            {
                double xyz[3];
                t8_forest_element_coordinate(forest, t, leaf, corner, xyz);
                const point p(xyz[0], xyz[1], xyz[2]);
                leafCorners[gc][corner] = p;
                c += p;
                getPoint(p); // register every corner now, so hanging-node
                             // midpoints (= finer cells' corners) are findable
            }
            c /= 8.0;
            leafCentre[gc] = c;
            // masterCell_ stores the OLD forest leaf index.  Convert it to the
            // corresponding OpenFOAM cell index after the first rebuild, where
            // leaf ordering and cell ordering are no longer guaranteed equal.
            label master = -1;
            if (gc < label(masterCell_.size()))
            {
                const label masterLeaf = masterCell_[gc];
                if
                (
                    masterLeaf >= 0
                 && masterLeaf < label(leafToCell_.size())
                )
                {
                    master = leafToCell_[masterLeaf];
                }
                else
                {
                    master = masterLeaf;
                }
            }
            meshMod.addCell(-1, -1, -1, master, -1);
        }
    }

    const double _tP1 = _ms();

    // ---- pass 2: faces -----------------------------------------------------
    for (t8_locidx_t t = 0; t < nTrees; ++t)
    {
        const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
        for (t8_locidx_t e = 0; e < nE; ++e)
        {
            const label gc = off[t] + e;
            const t8_element_t* leaf =
                t8_forest_get_leaf_element_in_tree(forest, t, e);
            const int myLevel =
                t8_element_get_level(scheme, T8_ECLASS_HEX, leaf);

            // 8 corner coords (z-order) -- cached from pass 1
            const std::array<point,8>& cc = leafCorners[gc];

            // 6 faces: corner indices + neighbour info
            int fck[6][4];
            int nn[6], nbl[6];
            label gnb[6];
            for (int k = 0; k < 6; ++k)
            {
                for (int j = 0; j < 4; ++j)
                    fck[k][j] = t8_element_get_face_corner
                        (scheme, T8_ECLASS_HEX, leaf, k, j);

                t8_element_t** neighs = nullptr;
                int* df = nullptr; int nN = 0;
                t8_locidx_t* ei = nullptr; t8_eclass_t ec;
                t8_forest_leaf_face_neighbors
                (
                    forest, t, leaf,
                    const_cast<const t8_element_t***>(&neighs),
                    k, &df, &nN, &ei, &ec
                );
                nn[k] = nN; gnb[k] = -1; nbl[k] = -1;
                if (nN >= 1)
                {
                    gnb[k] = ei[0];
                    nbl[k] = t8_element_get_level(scheme, ec, neighs[0]);
                }
                if (nN > 0) { T8_FREE(neighs); T8_FREE(ei); T8_FREE(df); }
            }

            for (int k = 0; k < 6; ++k)
            {
                bool emit = false, boundary = false;
                label own = gc, nei = -1;

                if (nn[k] == 0) { emit = true; boundary = true; }
                else if (nn[k] == 1)
                {
                    if (gnb[k] < nLeaf
                     && (nbl[k] < myLevel || (nbl[k]==myLevel && gc < gnb[k])))
                    { emit = true; own = min(gc,gnb[k]); nei = max(gc,gnb[k]); }
                }
                // nn[k] == 4: coarse side of a 2:1 face -> the finer cells emit

                if (!emit) continue;

                // polygon = the 4 face corners + any hanging-node midpoint that
                // exists on a face edge. A coarse edge's midpoint is in the
                // point set iff some (possibly diagonal) finer neighbour created
                // it as a corner -- this catches ALL 3D hanging nodes.
                std::vector<point> poly;
                for (int j = 0; j < 4; ++j) poly.push_back(cc[fck[k][j]]);
                for (int j1 = 0; j1 < 4; ++j1)
                {
                    for (int j2 = j1+1; j2 < 4; ++j2)
                    {
                        const int a = fck[k][j1], b = fck[k][j2];
                        if (__builtin_popcount(unsigned(a^b)) != 1) continue; // edge?
                        const point mid = 0.5*(cc[a]+cc[b]);
                        if (pointExists(mid)) poly.push_back(mid);
                    }
                }

                point fcen(0,0,0);
                for (const point& p : poly) fcen += p;
                fcen /= scalar(poly.size());
                const point ownC = leafCentre[own];
                const point oref = boundary ? ownC : leafCentre[nei];
                const label patch = boundary ? assignSide(fcen) : -1;
                addPolyFace(poly, own, nei, patch, ownC, oref, boundary);
            }
        }
    }

    const double _tP2 = _ms();

    // ---- commit the new mesh and map the fields ---------------------------
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);
    updateMesh(map);

    // cells were added in forest-leaf order (off[t]+e) and changeMesh did not
    // reorder them -> OF cell i == forest leaf i. Establish the identity map
    // used by the incremental rebuild.
    rebuildLeafToCell(); // robust forest-leaf <-> OF cell (centroid hash)

    const double _tCM = _ms();

    // patch sizes changed: re-evaluate BCs of all vol fields
    {
        HashTable<volScalarField*> sFlds =
            objectRegistry::lookupClass<volScalarField>();
        forAllIter(HashTable<volScalarField*>, sFlds, it)
        {
            const word& name = it.key();
            if (name != "nut")
            {
                (*it)->correctBoundaryConditions();
            }
        }
        HashTable<volVectorField*> vFlds =
            objectRegistry::lookupClass<volVectorField>();
        forAllIter(HashTable<volVectorField*>, vFlds, it)
            (*it)->correctBoundaryConditions();
    }

    // Keep the mapped (MULES-limited) fluxes; recompute from the
    // mapped velocity only faces without an old master (native rule).
    correctMappedFluxes(map);

    const double _tFin = _ms();

    Info<< "dynamicT8codeFvMesh: rebuilt OF mesh with " << nCells()
        << " cells, " << nFaces() << " faces, " << nPoints()
        << " points" << endl;
    Info<< "  rebuild breakdown [ms]: pass1(pts/cells)=" << _tP1
        << " pass2(faces)=" << _tP2 << " changeMesh+updateMesh=" << _tCM
        << " BC+flux=" << _tFin << endl;

    return true;
}


// Re-derive leafToCell_ (forest leaf -> OF cell) by locating each t8 leaf
// centroid in the rebuilt OF mesh. A quantised-centroid hash is not used
// here: coordinates that agree to 1e-16 can still round into different bins
// and silently scramble every mapped field on the next adapt.
void dynamicT8codeFvMesh::rebuildLeafToCell()
{
    t8_forest_t forest = FOREST;
    const t8_locidx_t nTrees = t8_forest_get_num_local_trees(forest);
    const t8_locidx_t nLeaf  = t8_forest_get_local_num_leaf_elements(forest);
    const labelList off = computeTreeOffsets(forest);

    leafToCell_.assign(nLeaf, -1);
    for (t8_locidx_t t = 0; t < nTrees; ++t)
    {
        const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
        for (t8_locidx_t e = 0; e < nE; ++e)
        {
            const t8_element_t* el =
                t8_forest_get_leaf_element_in_tree(forest, t, e);
            double xyz[3];
            t8_forest_element_centroid(forest, t, el, xyz);
            const point c
            (
                xyz[0], xyz[1],
                planar2D_ ? 0.5*(zMin_ + zMax_) : xyz[2]
            );
            const label cellI = this->findCell(c);
            if (cellI >= 0 && cellI < nCells())
            {
                leafToCell_[off[t] + e] = cellI;
            }
        }
    }
}


// Incremental rebuild: only the cells t8 changed (+ a finer-neighbour halo so
// the rebuilt region's boundary is clean) are removed and re-added; the rest of
// the OF mesh is retained. Reuses the from-scratch face logic on the region.
bool dynamicT8codeFvMesh::forestToOpenFOAMIncremental()
{
    t8_forest_t forest = FOREST;
    const t8_scheme_c* scheme = SCHEME;
    const t8_locidx_t nTrees = t8_forest_get_num_local_trees(forest);
    const t8_locidx_t nLeaf  = t8_forest_get_local_num_leaf_elements(forest);
    const labelList off = computeTreeOffsets(forest);

    auto _pt = std::chrono::steady_clock::now();
    auto _ms = [&]()
    {
        auto n = std::chrono::steady_clock::now();
        double d = std::chrono::duration<double,std::milli>(n - _pt).count();
        _pt = n; return d;
    };

    // ---- per-leaf element pointer / tree / level --------------------------
    std::vector<t8_locidx_t> leafTree(nLeaf);
    std::vector<const t8_element_t*> leafPtr(nLeaf);
    std::vector<int> leafLevel(nLeaf);
    for (t8_locidx_t t = 0; t < nTrees; ++t)
    {
        const t8_locidx_t nE = t8_forest_get_tree_num_leaf_elements(forest, t);
        for (t8_locidx_t e = 0; e < nE; ++e)
        {
            const label gc = off[t] + e;
            leafTree[gc] = t;
            leafPtr[gc]  = t8_forest_get_leaf_element_in_tree(forest, t, e);
            leafLevel[gc] = t8_element_get_level(scheme, T8_ECLASS_HEX, leafPtr[gc]);
        }
    }

    // neighbours of a leaf across face k: fills out[] with new-leaf indices,
    // returns count (0,1 or 4). lvlOut[] = neighbour levels.
    auto faceNbrs = [&](label gc, int k, label out[4], int lvlOut[4]) -> int
    {
        t8_element_t** neighs = nullptr; int* df = nullptr; int nN = 0;
        t8_locidx_t* ei = nullptr; t8_eclass_t ec;
        t8_forest_leaf_face_neighbors
        (
            forest, leafTree[gc], leafPtr[gc],
            const_cast<const t8_element_t***>(&neighs),
            k, &df, &nN, &ei, &ec
        );
        const int m = (nN > 4 ? 4 : nN);
        for (int i = 0; i < m; ++i)
        {
            out[i] = ei[i];
            lvlOut[i] = t8_element_get_level(scheme, ec, neighs[i]);
        }
        if (nN > 0) { T8_FREE(neighs); T8_FREE(ei); T8_FREE(df); }
        return nN;
    };

    // Cached neighbour query: each (leaf, face) hits t8_forest_leaf_face_neighbors
    // (which mallocs) at most ONCE -- the region BFS, the finer-closure and the
    // face pass all share this, instead of querying 2-3x. cNN<0 => not cached.
    std::vector<signed char> cNN(std::size_t(nLeaf)*6, -2);
    std::vector<label> cNb(std::size_t(nLeaf)*6*4);
    std::vector<signed char> cLv(std::size_t(nLeaf)*6*4);
    auto faceNbrsC = [&](label gc, int k, label out[4], int lvlOut[4]) -> int
    {
        const std::size_t b = (std::size_t(gc)*6 + k);
        if (cNN[b] == -2)
        {
            const int nN = faceNbrs(gc, k, out, lvlOut);
            cNN[b] = static_cast<signed char>(nN > 4 ? 4 : nN);
            const int m = (nN > 4 ? 4 : nN);
            for (int i = 0; i < m; ++i)
            { cNb[b*4+i] = out[i]; cLv[b*4+i] = static_cast<signed char>(lvlOut[i]); }
            return nN;
        }
        const int nN = cNN[b];
        for (int i = 0; i < nN; ++i)
        { out[i] = cNb[b*4+i]; lvlOut[i] = cLv[b*4+i]; }
        return nN;
    };

    // ---- Phase A: build the rebuild region --------------------------------
    // dirty cells, expanded by HALO=2 layers of ALL face-neighbours so that
    // every cell sharing a face OR edge with a level-changed cell is rebuilt
    // (a 3D edge is shared by 4 cells; the diagonal one is 2 face-hops away) ->
    // every NEW hanging node is emitted consistently on all faces touching it.
    // Then closed under "finer neighbour" so no region COARSE cell faces a kept
    // FINER cell (which it could not emit sub-faces for on its own).
    std::vector<unsigned char> inRegion(nLeaf, 0);
    std::vector<label> cur, nxt;
    for (t8_locidx_t i = 0; i < nLeaf; ++i)
        if (leafDirty_[i]) { inRegion[i] = 1; cur.push_back(i); }
    const int HALO = 2;
    for (int layer = 0; layer < HALO; ++layer)
    {
        nxt.clear();
        for (const label gc : cur)
            for (int k = 0; k < 6; ++k)
            {
                label nb[4]; int lv[4];
                const int nN = faceNbrsC(gc, k, nb, lv);
                for (int i = 0; i < nN; ++i)
                    if (nb[i] < nLeaf && !inRegion[nb[i]])
                    { inRegion[nb[i]] = 1; nxt.push_back(nb[i]); }
            }
        cur.swap(nxt);
    }
    std::vector<label> stack;
    for (t8_locidx_t i = 0; i < nLeaf; ++i) if (inRegion[i]) stack.push_back(i);
    while (!stack.empty())
    {
        const label gc = stack.back(); stack.pop_back();
        for (int k = 0; k < 6; ++k)
        {
            label nb[4]; int lv[4];
            const int nN = faceNbrsC(gc, k, nb, lv);
            for (int i = 0; i < nN; ++i)
                if (nb[i] < nLeaf && lv[i] > leafLevel[gc] && !inRegion[nb[i]])
                { inRegion[nb[i]] = 1; stack.push_back(nb[i]); }
        }
    }


    // ---- Phase B: removed OLD cells + kept new-cell ids -------------------
    std::vector<label> newCellOf(nLeaf, -1);
    std::vector<unsigned char> removedOld(leafToCell_.size(), 0);
    for (std::size_t L = 0; L < oldDirty_.size(); ++L)
        if (oldDirty_[L]) removedOld[leafToCell_[L]] = 1;
    for (t8_locidx_t gn = 0; gn < nLeaf; ++gn)
    {
        if (inRegion[gn])
        {
            if (!leafDirty_[gn])               // halo (clean) cell: rebuild it
                removedOld[leafToCell_[masterCell_[gn]]] = 1;
        }
        else                                    // kept: OF cell stays
            newCellOf[gn] = leafToCell_[masterCell_[gn]];
    }

    // ---- topology builder seeded with the CURRENT mesh -------------------
    directTopoChange meshMod(*this);

    for (label c = 0; c < nCells(); ++c)
        if (removedOld[c]) meshMod.removeCell(c, -1);
    {
        const labelUList& fOwn = faceOwner();
        const labelUList& fNei = faceNeighbour();
        const label nInt = nInternalFaces();
        forAll(faces(), f)
        {
            const label own = fOwn[f];
            const label nei = (f < nInt ? fNei[f] : -1);
            if (removedOld[own] || (nei >= 0 && removedOld[nei]))
                meshMod.removeFace(f, -1);
        }
    }

    // ---- point dedup: pre-load with ALL current mesh points --------------
    std::unordered_map<std::array<long,3>, label, Key3Hash> pointMap;
    const double qtol = 1e-9;
    auto pkey = [&](const point& p)
    {
        return std::array<long,3>
            { std::lround(p.x()/qtol), std::lround(p.y()/qtol),
              std::lround(p.z()/qtol) };
    };
    {
        const pointField& mp = points();
        pointMap.reserve(mp.size()*2);
        forAll(mp, pI) pointMap[pkey(mp[pI])] = pI;
    }
    const label origNPts = points().size();
    std::vector<char> regCorner(origNPts, 0); // existing point is a region corner
    auto getPoint = [&](const point& p) -> label
    {
        const auto k = pkey(p);
        auto it = pointMap.find(k);
        if (it != pointMap.end()) return it->second;
        const label id = meshMod.addPoint(p, -1, -1, true);
        pointMap[k] = id;
        return id;
    };
    auto pointExists = [&](const point& p) -> bool
    { return pointMap.find(pkey(p)) != pointMap.end(); };

    // ---- domain-side patch assignment (6 sides) --------------------------
    const boundBox& bb = bounds();
    const point bMin = bb.min(), bMax = bb.max();
    auto sideOf = [&](const point& p) -> int
    {
        const scalar d[6] = { mag(p.x()-bMin.x()), mag(p.x()-bMax.x()),
            mag(p.y()-bMin.y()), mag(p.y()-bMax.y()),
            mag(p.z()-bMin.z()), mag(p.z()-bMax.z()) };
        int s = 0; for (int i = 1; i < 6; ++i) if (d[i] < d[s]) s = i; return s;
    };
    const polyBoundaryMesh& pbm = boundaryMesh();
    labelList sidePatchID(6, -1);
    forAll(pbm, pI)
    {
        const vectorField& fc = pbm[pI].faceCentres();
        point c(0,0,0); forAll(fc, k) c += fc[k];
        if (fc.size()) c /= fc.size();
        sidePatchID[sideOf(c)] = pI;
    }

    // OLD-mesh face data (valid until changeMesh below) for face masters:
    // each new face maps from the old face with the closest centre among
    // its cells' old masters, so mapped surface fields stay finite (see
    // the from-scratch branches: unmapped faces left e.g. nHatf_
    // uninitialized -> NaN -> SIGFPE in the flux weights).
    const cellList& oldCellFaces = cells();
    vectorField oldFC(nFaces());
    {
        const vectorField& ifc = Cf().internalField();
        forAll(ifc, f) oldFC[f] = ifc[f];
        const polyBoundaryMesh& obm = boundaryMesh();
        forAll(obm, p)
        {
            const vectorField& pfc = obm[p].faceCentres();
            const label s = obm[p].start();
            forAll(pfc, i) oldFC[s + i] = pfc[i];
        }
    }
    // New-cell index -> old-cell master (filled in Phase C below). Kept
    // cells are absent: their pre-changeMesh index already is an old cell.
    std::unordered_map<label, label> newCellMaster;

    auto addPolyFace = [&](std::vector<point> poly, label own, label nei,
                           label patch, const point& ownC, const point& oref,
                           bool boundary)
    {
        const int np = int(poly.size());
        if (np < 3) return;
        point mn = poly[0], mx = poly[0];
        for (const point& p : poly) { mn = min(mn,p); mx = max(mx,p); }
        const vector ext = mx - mn;
        int axis = 0; if (ext[1] < ext[axis]) axis = 1; if (ext[2] < ext[axis]) axis = 2;
        const int a1 = (axis+1)%3, a2 = (axis+2)%3;
        point ctr(0,0,0); for (const point& p : poly) ctr += p; ctr /= scalar(np);
        std::sort(poly.begin(), poly.end(), [&](const point& A, const point& B)
            { return std::atan2(A[a2]-ctr[a2],A[a1]-ctr[a1])
                   < std::atan2(B[a2]-ctr[a2],B[a1]-ctr[a1]); });
        vector nrm(0,0,0);
        for (int i = 0; i < np; ++i)
        {
            const point& P = poly[i]; const point& Q = poly[(i+1)%np];
            nrm.x() += (P.y()-Q.y())*(P.z()+Q.z());
            nrm.y() += (P.z()-Q.z())*(P.x()+Q.x());
            nrm.z() += (P.x()-Q.x())*(P.y()+Q.y());
        }
        const vector dir = boundary ? (ctr-ownC) : (oref-ownC);
        const bool rev = ((nrm & dir) < 0);
        // Face master: old face with the closest centre among the two
        // cells' old masters (kept cells map to themselves).
        label masterFace = -1;
        {
            const label nwc[2] = {own, nei};
            label mc[2] = {-1, -1};
            for (int k = 0; k < 2; ++k)
            {
                const label nc = nwc[k];
                if (nc < 0) continue;
                label oc = -1;
                auto it = newCellMaster.find(nc);
                if (it != newCellMaster.end())
                {
                    oc = it->second;
                }
                else
                {
                    oc = nc;
                }
                if (oc >= 0 && oc < oldCellFaces.size())
                {
                    mc[k] = oc;
                }
            }
            scalar bestD = GREAT;
            for (int k = 0; k < 2; ++k)
            {
                if (mc[k] < 0) continue;
                if (k == 1 && mc[1] == mc[0]) continue;
                const cell& cFaces = oldCellFaces[mc[k]];
                forAll(cFaces, fi)
                {
                    const label of = cFaces[fi];
                    if (of < 0 || of >= oldFC.size()) continue;
                    const scalar d = mag(oldFC[of] - ctr);
                    if (d < bestD)
                    {
                        bestD = d;
                        masterFace = of;
                    }
                }
            }
        }
        labelList fp(np);
        for (int i = 0; i < np; ++i) fp[i] = getPoint(rev ? poly[np-1-i] : poly[i]);
        meshMod.addFace
        (
            face(fp), own, nei, -1, -1, masterFace,
            false, patch, -1, false
        );
    };

    // ---- Phase C: add region cells (+corners) ----------------------------
    std::vector<std::array<point,8>> rc(nLeaf);  // corners of region leaves
    for (t8_locidx_t gn = 0; gn < nLeaf; ++gn)
    {
        if (!inRegion[gn]) continue;
        const label t = leafTree[gn]; const t8_element_t* leaf = leafPtr[gn];
        point ctr(0,0,0);
        for (int corner = 0; corner < 8; ++corner)
        {
            double xyz[3];
            t8_forest_element_coordinate(forest, t, leaf, corner, xyz);
            const point p(xyz[0],xyz[1],xyz[2]);
            const label pid = getPoint(p);
            if (pid < origNPts) regCorner[pid] = 1;
            rc[gn][corner] = p; ctr += p;
        }
        (void)ctr;
        const label master = leafToCell_[masterCell_[gn]];
        newCellOf[gn] = meshMod.addCell(-1,-1,-1, master, -1);
        newCellMaster[newCellOf[gn]] = master;
    }

    auto leafCtr = [&](label gn) -> point
    {
        if (inRegion[gn])
        { point c(0,0,0); for (int i=0;i<8;i++) c += rc[gn][i]; return c/8.0; }
        // kept cell centre
        return this->C().internalField()[ newCellOf[gn] ];
    };

    // ---- Phase D: faces of the region ------------------------------------
    const double _tPrep = _ms();
    for (t8_locidx_t gn = 0; gn < nLeaf; ++gn)
    {
        if (!inRegion[gn]) continue;
        const std::array<point,8>& cc = rc[gn];
        int fckArr[6][4];
        for (int k = 0; k < 6; ++k)
            for (int j = 0; j < 4; ++j)
                fckArr[k][j] = t8_element_get_face_corner(scheme, T8_ECLASS_HEX, leafPtr[gn], k, j);

        for (int k = 0; k < 6; ++k)
        {
            label nb[4]; int lv[4];
            const int nN = faceNbrsC(gn, k, nb, lv);

            bool emit = false, boundary = false;
            label gNbr = -1;                 // neighbour leaf (-1 => boundary)

            if (nN == 0) { emit = true; boundary = true; }
            else if (nN == 1)
            {
                const label g1 = nb[0];
                if (g1 < nLeaf)
                {
                    if (!inRegion[g1])              // neighbour kept -> always emit
                    { emit = true; gNbr = g1; }
                    else if (lv[0] < leafLevel[gn]) // I'm finer -> emit
                    { emit = true; gNbr = g1; }
                    else if (lv[0] == leafLevel[gn] && gn < g1)
                    { emit = true; gNbr = g1; }
                }
            }
            // nN == 4: I'm the coarse side; all finer neighbours are in the
            // region (Phase A) and emit their own sub-faces -> skip.

            if (!emit) continue;

            std::vector<point> poly;
            for (int j = 0; j < 4; ++j) poly.push_back(cc[fckArr[k][j]]);
            for (int j1 = 0; j1 < 4; ++j1)
              for (int j2 = j1+1; j2 < 4; ++j2)
              {
                  const int a = fckArr[k][j1], b = fckArr[k][j2];
                  if (__builtin_popcount(unsigned(a^b)) != 1) continue;
                  const point mid = 0.5*(cc[a]+cc[b]);
                  if (pointExists(mid)) poly.push_back(mid);
              }
            point fcen(0,0,0); for (const point& p : poly) fcen += p;
            fcen /= scalar(poly.size());

            // OpenFOAM requires owner < neighbour for internal faces. Order the
            // two cells and orient the face from owner to neighbour.
            label own, nei; point ownC, oref;
            const label patch = boundary ? sidePatchID[sideOf(fcen)] : -1;
            if (boundary)
            {
                if (patch < 0) continue;
                own = newCellOf[gn]; nei = -1; ownC = leafCtr(gn); oref = ownC;
            }
            else
            {
                const label cThis = newCellOf[gn], cNbr = newCellOf[gNbr];
                const point cThisC = leafCtr(gn), cNbrC = leafCtr(gNbr);
                if (cThis < cNbr)
                { own = cThis; nei = cNbr; ownC = cThisC; oref = cNbrC; }
                else
                { own = cNbr; nei = cThis; ownC = cNbrC; oref = cThisC; }
            }
            addPolyFace(poly, own, nei, patch, ownC, oref, boundary);
        }
    }
    // ---- remove points orphaned by the removed (coarsened) region --------
    // A point disappears iff every cell using it was removed AND it is not a
    // corner of any re-added region cell.
    {
        const labelListList& pc = pointCells();
        for (label p = 0; p < origNPts; ++p)
        {
            if (regCorner[p]) continue;
            const labelList& cs = pc[p];
            bool allRemoved = true;
            forAll(cs, i)
                if (!removedOld[cs[i]]) { allRemoved = false; break; }
            if (allRemoved) meshMod.removePoint(p, -1);
        }
    }

    const double _tFaces = _ms();

    // ---- commit + map fields + refresh BCs/fluxes ------------------------
    // orderCells=true: bandCompress the cell numbering. The retained+appended
    // ordering of an incremental change has poor locality, which degrades the
    // linear solvers; a from-scratch mesh is naturally SFC-ordered.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false, true, true);
    updateMesh(map);
    const double _tCM = _ms();

    {
        HashTable<volScalarField*> sFlds = objectRegistry::lookupClass<volScalarField>();
        forAllIter(HashTable<volScalarField*>, sFlds, it)
        {
            const word& name = it.key();
            if (name != "nut")
            {
                (*it)->correctBoundaryConditions();
            }
        }
        HashTable<volVectorField*> vFlds = objectRegistry::lookupClass<volVectorField>();
        forAllIter(HashTable<volVectorField*>, vFlds, it) (*it)->correctBoundaryConditions();
    }

    // Keep the mapped (MULES-limited) fluxes; recompute from the
    // mapped velocity only faces without an old master (native rule).
    correctMappedFluxes(map);

    // re-establish forest leaf <-> OF cell for the next adapt
    rebuildLeafToCell();

    label nRegion = 0; for (t8_locidx_t i=0;i<nLeaf;i++) nRegion += inRegion[i];
    Info<< "dynamicT8codeFvMesh: incremental rebuild, region=" << nRegion
        << "/" << label(nLeaf) << " cells, now " << nCells() << " cells, "
        << nFaces() << " faces, " << nPoints() << " points" << nl
        << "  incr breakdown [ms]: prep=" << _tPrep << " faces=" << _tFaces
        << " changeMesh=" << _tCM << endl;
    return true;
}


void dynamicT8codeFvMesh::copyFieldsFromForest()
{
    // Write the transferred per-leaf values onto the rebuilt OpenFOAM fields.
    // changeMesh() can renumber cells, so never assume OF cell i == t8 leaf i;
    // rebuildLeafToCell() establishes the required leaf -> cell map.
    const label nLeaf = alpha1Data_.size();
    const bool useMap = (label(leafToCell_.size()) == nLeaf);
    auto cellOf = [&](label leaf) -> label
    { return useMap ? leafToCell_[leaf] : leaf; };

    if (foundObject<volScalarField>(fieldName_))
    {
        volScalarField& a =
            const_cast<volScalarField&>(lookupObject<volScalarField>(fieldName_));
        for (label leaf = 0; leaf < nLeaf; ++leaf)
        {
            const label c = cellOf(leaf);
            if (c >= 0 && c < a.size()) a[c] = alpha1Data_[leaf];
        }
        a.correctBoundaryConditions();
    }
    if (foundObject<volVectorField>("U") && label(uData_.size()) >= 3*nLeaf)
    {
        volVectorField& U =
            const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
        for (label leaf = 0; leaf < nLeaf; ++leaf)
        {
            const label c = cellOf(leaf);
            if (c >= 0 && c < U.size())
            {
                U[c].x() = uData_[3*leaf + 0];
                U[c].y() = uData_[3*leaf + 1];
                U[c].z() = uData_[3*leaf + 2];
            }
        }
        U.correctBoundaryConditions();
    }
    if (foundObject<volScalarField>("pd") && label(pdData_.size()) >= nLeaf)
    {
        volScalarField& pd =
            const_cast<volScalarField&>(lookupObject<volScalarField>("pd"));
        for (label leaf = 0; leaf < nLeaf; ++leaf)
        {
            const label c = cellOf(leaf);
            if (c >= 0 && c < pd.size()) pd[c] = pdData_[leaf];
        }
        pd.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool dynamicT8codeFvMesh::update()
{
    // Run at most once per time step
    if (curTimeIndex_ == time().timeIndex())
    {
        return false;
    }
    curTimeIndex_ = time().timeIndex();

    if (!forestInitialized_)
    {
        initForest();
        return false;
    }

    if
    (
        refineInterval_ <= 0
     || time().timeIndex() == 0
     || (time().timeIndex() % refineInterval_) != 0
    )
    {
        return false;
    }

    using clk = std::chrono::steady_clock;
    const auto tStart = clk::now();

    const bool changed = runT8codeAdapt();
    const auto tAdapt = clk::now();

    bool meshChanged = false;
    if (changed)
    {
        // First rebuild (or non-incremental mode) is from-scratch and also
        // establishes leafToCell_; subsequent ones can be incremental.
        if (incrementalRebuild_ && !leafToCell_.empty())
        {
            meshChanged = forestToOpenFOAMIncremental();
        }
        else
        {
            meshChanged = forestToOpenFOAM();
        }

        if (meshChanged)
        {
            copyFieldsFromForest();
        }
    }
    const auto tEnd = clk::now();

    auto ms = [](clk::time_point a, clk::time_point b)
    {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    Info<< "dynamicT8codeFvMesh::update timing [ms]: t8 adapt = "
        << ms(tStart, tAdapt) << ", OF rebuild = " << ms(tAdapt, tEnd)
        << ", total = " << ms(tStart, tEnd) << endl;

    return meshChanged;
}


// ************************************************************************* //
