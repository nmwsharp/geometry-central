#pragma once

#include "geometrycentral/geometry.h"
#include "geometrycentral/halfedge_mesh.h"

#include <deque>
#include <vector>

namespace geometrycentral {

// Used below
struct SegmentEndpoint {

  // Construct an edge crossing endpoint
  SegmentEndpoint(HalfedgePtr he, double tCross_) {
    isEdgeCrossing = true;
    halfedge = he;
    tCross = tCross_;
  }

  // Construct a face endpoint
  SegmentEndpoint(FacePtr face_, Vector3 faceCoords_) {
    isEdgeCrossing = false;
    face = face_;
    faceCoords = faceCoords_;
  }

  // Most of the time, this is an edge crossing
  bool isEdgeCrossing;
  HalfedgePtr halfedge; // enters at halfedges, exits halfedge.twin()
  double tCross;        // along halfedge

  // Can also be a point in a face
  FacePtr face;
  Vector3 faceCoords;
};

// A piecewise linear chunk within a face
struct CurveSegment {
  FacePtr face;

  Vector3 startBaryCoord;
  Vector3 endBaryCoord;

  Vector3 startPosition;
  Vector3 endPosition;

  HalfedgePtr startHe; // both inside of the face
  HalfedgePtr endHe;


  double length();
};


// A curve embedded in the surface of a mesh. Defined to be linear within every face, stored as a sequence of edge
// crossings with possible endpoints in faces. May be closed or open. May or may not have self-crossings.
// The representation is oriented, but mostly you can ignore this to work with unorient curves.
class MeshEmbeddedCurve {

public:
  MeshEmbeddedCurve(Geometry<Euclidean>* geometry);
  MeshEmbeddedCurve(){};

  // Build the curve
  // void extendFront(FacePtr f, Vector3 bCoord);
  void extendBack(FacePtr f, Vector3 bCoord);
  void tryExtendBack(FacePtr f,
                     Vector3 bCoord); // same as extend back, but no-ops if face is not adjacent to current end of curve
  void removeLastEndpoint();
  void removeFirstEndpoint();
  void rotateArbitraryStart(); // for a closed curve, shift the "seam" forward along the curve
  void closeCurve();
  void clearCurve();

  void setFromZeroLevelset(VertexData<double>& implicitF);


  // Get data about the curve
  bool isClosed();
  std::vector<CurveSegment> getCurveSegments();
  FacePtr endingFace(bool reportForClosed=false);
  FacePtr startingFace(bool reportForClosed=false);
  double computeLength();
  size_t nSegments();
  
  bool crossesFace(FacePtr f);

  // Throws an error if this is not a valid closed or open curve. Does not check self-intersection.
  void validate();

  // Copy this curve on to another equivalent mesh object using a transfer map. Here the transfer map is defined to go
  // from otherGeom.getMesh() to this::mesh
  MeshEmbeddedCurve copy(HalfedgeMeshDataTransfer& transfer, Geometry<Euclidean>* otherGeom);
  // Copies in teh opposite direction
  MeshEmbeddedCurve copyBack(HalfedgeMeshDataTransfer& transfer, Geometry<Euclidean>* otherGeom);

private:
  HalfedgeMesh* mesh = nullptr;
  Geometry<Euclidean>* geometry = nullptr;

  std::deque<SegmentEndpoint> segmentPoints;

  // == Helpers
  Vector3 barycoordsForHalfedgePoint(HalfedgePtr he, double t);
  Vector3 positionOfSegmentEndpoint(SegmentEndpoint& p);
  FacePtr faceBefore(SegmentEndpoint& p);
  FacePtr faceAfter(SegmentEndpoint& p);
  HalfedgePtr connectingHalfedge(FacePtr f1, FacePtr f2);
  bool facesAreAdjacentOrEqual(FacePtr f1, FacePtr f2);
  double crossingPointAlongEdge(HalfedgePtr sharedHe, Vector3 bCoord1, Vector3 bCoord2);
  double scalarFunctionZeroPoint(double f0, double f1);
};
}
