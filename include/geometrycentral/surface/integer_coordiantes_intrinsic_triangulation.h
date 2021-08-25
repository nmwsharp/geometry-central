#pragma once

#include "geometrycentral/surface/barycentric_coordinate_helpers.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/intrinsic_mollification.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <deque>

#include "CommonSubdivision.h"
#include "NormalCoordinates.h"

#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

struct CompoundCurve {
    std::vector<Curve> components;
};
inline std::ostream& operator<<(std::ostream& out, const CompoundCurve& curve) {
    out << "CompoundCurve {components: " << curve.components << "}";
    return out;
}

class NormalCoordinateIntrinsicTriangulation
    : public IntrinsicGeometryInterface {

  public:
    // Construct an intrinsic triangulation which sits atop this input mesh.
    // Initially, the input triangulation will just be a copy of the input mesh.
    NormalCoordinateIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom);

    // ======================================================
    //                   Core Members
    // ======================================================
    // The underlying surface on which the intrinsic triangulation has been
    // constructed
    ManifoldSurfaceMesh& inputMesh;
    IntrinsicGeometryInterface& inputGeom;

    // The connectivity of the intrinsic triangulation
    // note that somewhat confusingly, there is a .mesh reference which points
    // to this same mesh, inherited from the geometry interface
    std::unique_ptr<ManifoldSurfaceMesh> intrinsicMesh;

    std::unique_ptr<CommonSubdivision> commonSubdivision;

    EdgeData<double> intrinsicEdgeLengths; // length of each edge

    // TODO: At this point, should this just be a member variable? Why do I even
    // use unique_ptrs?
    std::unique_ptr<NormalCoordinates> normalCoordinates_ptr;
    NormalCoordinates& normalCoordinates;

    // location of vertices of intrinsicMesh as SurfacePoints on inputMesh
    VertexData<SurfacePoint> intrinsicToInput, inputToIntrinsic;

    // TODO: move outside class? - maybe delaunay flips should be external too
    // Edges that are not allowed to be flipped. In delaunay refinement, these
    // are split instead
    // std::vector<Edge> fixedEdges; // NSHARP: removed this since it was unused
    // and annoying to maintain
    EdgeData<bool> isFixed;

    double tracingTime     = 0;
    double insertionTime   = 0;
    size_t easyInsertions1 = 0;
    size_t easyInsertions2 = 0;
    size_t hardInsertions  = 0;
    size_t faceSplits      = 0;
    size_t edgeSplits      = 0;
    FaceData<bool> hasBeenSplit;

    // ======================================================
    //                 Queries & Accessors
    // ======================================================

    // Trace out the edges of the intrinsic triangulation along the surface of
    // the input mesh. Each path is ordered along edge.halfedge(), and includes
    // both the start and end points
    CommonSubdivision& traceEdges(bool geodesic = true);

    // Return by reference to make tests easier to write
    void traceEdges(CommonSubdivision& cs, bool geodesic = true);

    std::vector<std::vector<SurfacePoint>> traceEdgeSet();

    // ======================================================
    //                High-Level Mutators
    // ======================================================

    // Flips edge in the intrinsic triangulation until is satisfies the
    // intrinsic Delaunay criterion
    // Returns the number of flips performed
    size_t flipToDelaunay(double tol = 1e-5);

    // ======================================================
    //                Low-Level Mutators
    // ======================================================
    //
    // Basic operations to modify the intrinsic triangulation
    // NOTE: The individual operations do not call refreshQuantities(), so you
    // should call it if you want quantities updated.

    // If the edge is not Delaunay, flip it. Returns true if flipped.
    bool flipEdgeIfNotDelaunay(Edge e, double possibleEPS = 1e-5);

    // If the edge can be flipped, flip it (must be combinatorially flippable
    // and inside a convex quad). Returns true if flipped.
    bool flipEdgeIfPossible(Edge e, double possibleEPS = 1e-5,
                            bool verbose = false);
    double checkFlip(Edge e);

    // Insert circumcenter or split segment, from NIT, geometrycentral
    // TODO: do edge encroachment properly?
    Vertex insertCircumcenterOrSplitSegment(Face f, bool verbose = false);

    Vertex insertVertex(SurfacePoint pt, bool verbose = false);
    Vertex splitFace(Face f, Vector3 bary, bool verbose = false);

    Vertex splitEdge(Edge e, double bary, bool verbose = false);
    Vertex splitInteriorEdge(Edge e, double bary, bool verbose = false);
    Vertex splitBoundaryEdge(Edge e, double bary, bool verbose = false);

    // Try to remove an inserted vertex
    // Returns Face() if anything goes wrong
    Face removeInsertedVertex(Vertex v);


    // Move a vertex `v` in direction `vec`, represented as a vector in the
    // vertex's tangent space.
    Vertex moveVertex(Vertex v, Vector2 vec);

    // Assumes intrinsicEdgeLengths is up to date
    void updateCornerAngle(Corner c);

    // Assumes cornerAngles exists and is up to date
    void updateVertexAngleSum(Vertex v);

    // Assumes that intrinsicEdgeLengths is up to date
    void updateFaceArea(Face f);

    // Assumes cornerAngles, vertexAngleSums exist and are up to date
    void updateHalfedgeVectorsInVertex(Vertex v);

    // Assumes intrinsicEdgeLengths, faceAreas exist and are up to date
    void updateHalfedgeVectorsInFace(Face f);

    // ======================================================
    //                Low-Level Queries
    // ======================================================
    // TODO: use geometry-central properly
    double getCornerAngle(Corner c) const;
    double getCornerAngleCosine(Corner c) const;
    double getFaceArea(Face f) const;
    double getHalfedgeCotanWeight(Halfedge he) const;
    double getEdgeCotanWeight(Edge e) const;

    double getCircumradius(Face f) const;
    double getShortestEdge(Face f) const;

    double getMinAngle() const;
    double getMaxCircumcircleToEdgeLengthRatio() const;
    double getMinInteriorAngleSum() const;
    double getMinBoundaryAngleSum() const;

    // Get min angle on faces whose vertices all have angle sum at least
    // minAngleSum
    double getMinAngleOnValidFaces(double minAngleSum) const;

    // Get max cimrcumcircle-to-edge-length ratio on faces whose vertices all
    // have angle sum at least minAngleSum
    double getMaxRatioOnValidFaces(double minAngleSum) const;

    // Check if all face vertices have angle sum at least minAngleSum
    bool faceHasLargeAngleSums(Face f, double minAngleSum) const;

    // If f is contained entirely in some input face, check that its vertices
    // have angle sum at least minAngleSum
    // Otherwise return true
    bool parentFaceHasLargeAngleSums(Face f, double minAngleSum) const;

    bool isOnFixedEdge(Vertex v) const;

    // Takes in a halfedge of the intrinsic mesh whose edge's normal coordinate
    // is negative (meaning that it lies along an edge of the input mesh) and
    // returns the halfedge in the input mesh pointing in the same direction
    // e.vertex() must live in both meshes
    Halfedge getSharedInputEdge(Halfedge e) const;

    bool isDelaunay(Edge e, double tol) const;
    // bool isDelaunay(Edge e, double flippedLength);

    // Compute the number of vertices in the common subdivision
    // i.e. intrinsicMesh->nVertices() plus the sum of the normal coordinates
    size_t nSubdividedVertices() const;

    // HACK: represents arcs parallel to a mesh edge with a single pair {-n,
    // he} where n is the number of arcs parallel to he.edge() Trace an edge
    // of the input mesh over the intrinsic triangulation
    CompoundCurve traceInputEdge(Edge e, bool verbose = false) const;

    std::pair<bool, Curve> traceNextCurve(const Curve& oldCurve,
                                          bool verbose = false) const;

    // Inverse to traceInputEdge
    Halfedge identifyInputEdge(const Curve& path, bool verbose = false) const;

    // Identify shared halfedge, throw exception if halfedge is not shared
    // (i.e. edgeCoords[he.edge()] must be negative)
    Halfedge identifyInputEdge(Halfedge he) const;

    std::array<Vector2, 3> vertexCoordinatesInFace(Face face) const;

    void setFixedEdges(const EdgeData<bool>& fixedEdges);

    // If f is entirely contained in some face of the input mesh, return that
    // face Otherwise return Face()
    Face getParentFace(Face f) const;

    // ======================================================
    //                     Callbacks
    // ======================================================
    //
    // Get called whenever mesh mutations occur. Register a callback by
    // inserting it in to this list.

    // edge E if flipped
    std::list<std::function<void(Edge)>> edgeFlipCallbackList;

    // old face F is split by new vertex V
    std::list<std::function<void(Face, Vertex)>> faceSplitCallbackList;

    // old edge E is split to halfedge HE1,HE2 both with he.vertex() as split
    // vertex
    std::list<std::function<void(Edge, Halfedge, Halfedge)>>
        edgeSplitCallbackList;

  private:
    // ======================================================
    //               Geometry Interface
    // ======================================================
    //
    // Satisfy the requirements of the IntrinsicGeometryInterface

    // Override the compute edge lengths method from intrinsic geometry.
    virtual void computeEdgeLengths() override;

    // Isometrically lay out the vertices around a halfedge in 2D
    // coordinates he points from vertex 2 to 0; others are numbered CCW
    std::array<Vector2, 4> layoutDiamond(Halfedge he) const;

    // Callback helpers
    void invokeEdgeFlipCallbacks(Edge e);
    void invokeFaceSplitCallbacks(Face f, Vertex v);
    void invokeEdgeSplitCallbacks(Edge e, Halfedge he1, Halfedge he2);
};

// Compute the cotan weight of halfedge ij in terms of the lengths of its
// neighbors
double halfedgeCotanWeight(double lij, double ljk, double lki);

FaceData<Vector2>
interpolateTangentVectorsB(const NormalCoordinateIntrinsicTriangulation& tri,
                           const CommonSubdivision& cs,
                           const FaceData<Vector2>& dataB);

#include "NormalCoordinateIntrinsicTriangulation.ipp"
