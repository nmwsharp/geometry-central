#include "geometrycentral/mesh_embedded_curve.h"

using namespace geometrycentral;

using std::cout;
using std::endl;

MeshEmbeddedCurve::MeshEmbeddedCurve(Geometry<Euclidean>* geometry_) : geometry(geometry_) {

  mesh = geometry->getMesh();
  if (!mesh->isSimplicial()) {
    throw std::runtime_error("Embedded curves only supported on simplicial mesh");
  }
}


// Helper utilities
Vector3 MeshEmbeddedCurve::barycoordsForHalfedgePoint(HalfedgePtr he, double t) {

  Vector3 bCoord{0.0, 0.0, 0.0};
  HalfedgePtr currHe = he.face().halfedge();
  for (size_t i = 0; i < 3; i++) {
    if (currHe == he) {
      bCoord[i] = 1.0 - t;
      bCoord[(i + 1) % 3] = t;
    }
    currHe = currHe.next();
  }

  return bCoord;
}


Vector3 MeshEmbeddedCurve::positionOfSegmentEndpoint(SegmentEndpoint& p) {

  if (p.isEdgeCrossing) {
    return geometry->position(p.halfedge.vertex()) + p.tCross * geometry->vector(p.halfedge);
  } else {
    Vector3 pos{0.0, 0.0, 0.0};
    HalfedgePtr currHe = p.face.halfedge();
    for (size_t i = 0; i < 3; i++) {
      pos += geometry->position(currHe.vertex()) * p.faceCoords[i];
      currHe = currHe.next();
    }
    return pos;
  }
}

FacePtr MeshEmbeddedCurve::faceBefore(SegmentEndpoint& p) {
  if (p.isEdgeCrossing) {
    return p.halfedge.face();
  } else {
    return p.face;
  }
}

FacePtr MeshEmbeddedCurve::faceAfter(SegmentEndpoint& p) {
  if (p.isEdgeCrossing) {
    return p.halfedge.twin().face();
  } else {
    return p.face;
  }
}

HalfedgePtr MeshEmbeddedCurve::connectingHalfedge(FacePtr f1, FacePtr f2) {
  for (HalfedgePtr he : f1.adjacentHalfedges()) {
    if (he.twin().face() == f2) {
      return he;
    }
  }
  throw std::runtime_error("Faces do not share an adjacent halfedge");
}

bool MeshEmbeddedCurve::facesAreAdjacentOrEqual(FacePtr f1, FacePtr f2) {
  if (f1 == f2) return true;
  for (HalfedgePtr he : f1.adjacentHalfedges()) {
    if (he.twin().face() == f2) {
      return true;
    }
  }
  return false;
}

double MeshEmbeddedCurve::crossingPointAlongEdge(HalfedgePtr sharedHe, Vector3 bCoord1, Vector3 bCoord2) {


  // Build a coordinate space with the shared edge as the y coordinate
  Vector3 basisY = geometry->vector(sharedHe);
  double bLen = norm(basisY);
  basisY /= bLen;
  Vector3 rootP = geometry->position(sharedHe.vertex());

  // Coordinates in first face
  Vector2 p1{0.0, 0.0};
  Vector3 basisX1 = cross(basisY, geometry->normal(sharedHe.face()));
  HalfedgePtr currHe = sharedHe.face().halfedge();
  for (size_t i = 0; i < 3; i++) {
    Vector3 pv = geometry->position(currHe.vertex()) - rootP;
    p1 += Vector2{dot(basisX1, pv), dot(basisY, pv)} * bCoord1[i];
    currHe = currHe.next();
  }
  p1 /= bLen;

  // Coordinates in second face
  Vector2 p2{0.0, 0.0};
  Vector3 basisX2 = cross(basisY, geometry->normal(sharedHe.twin().face()));
  currHe = sharedHe.twin().face().halfedge();
  for (size_t i = 0; i < 3; i++) {
    Vector3 pv = geometry->position(currHe.vertex()) - rootP;
    p2 += Vector2{dot(basisX2, pv), dot(basisY, pv)} * bCoord2[i];
    currHe = currHe.next();
  }
  p2 /= bLen;

  // lol slope intercept
  double slope = (p2.y - p1.y) / (p2.x - p1.x);
  double intercept = (p2.y - slope * p2.x);

  if (intercept < -1e-3 || intercept > 1 + 1e-3) throw std::runtime_error("crossing calculation failed");

  return intercept;
}

double MeshEmbeddedCurve::scalarFunctionZeroPoint(double f0, double f1) {
  double EPS = 1e-4;
  return clamp(f0 / (f0 - f1), EPS, 1.0 - EPS);
}


void MeshEmbeddedCurve::tryExtendBack(FacePtr f, Vector3 bCoord) {
  if (segmentPoints.size() == 0 || facesAreAdjacentOrEqual(f, faceAfter(segmentPoints.back()))) {
    extendBack(f, bCoord);
  }
}

void MeshEmbeddedCurve::extendBack(FacePtr f, Vector3 bCoord) {

  // Special case if this is the start of the curve
  if (segmentPoints.size() == 0) {
    segmentPoints.push_back(SegmentEndpoint(f, bCoord));
    return;
  }

  if (isClosed()) {
    throw std::runtime_error("Can't extend a closed curve");
  }

  SegmentEndpoint oldEnd = segmentPoints.back();

  // If the new point is the same face as the old, just move it
  if (oldEnd.face == f) {
    segmentPoints.back().faceCoords = bCoord;
  }
  // The new point is in an adjacent face
  else {

    // Handle the case where we're just starting to construct a curve, and the previous point is a first face point
    if (segmentPoints.size() > 1) {
      segmentPoints.pop_back();
    }

    // Create a new crossing
    HalfedgePtr sharedHe = connectingHalfedge(faceAfter(oldEnd), f);
    double t = crossingPointAlongEdge(sharedHe, oldEnd.faceCoords, bCoord);
    segmentPoints.push_back(SegmentEndpoint(sharedHe, t));

    // Create a new endpoint
    segmentPoints.push_back(SegmentEndpoint(f, bCoord));
  }
}

void MeshEmbeddedCurve::removeFirstEndpoint() {

  if (segmentPoints.size() == 0) return;

  if (segmentPoints.size() == 1) {
    segmentPoints.pop_back();
    return;
  }

  // Break up a closed curve
  if (isClosed()) {

    // Remove the old crossing
    SegmentEndpoint crossPoint = segmentPoints.front();
    segmentPoints.pop_front();

    // Create a new face point in the end face
    segmentPoints.push_back(SegmentEndpoint(crossPoint.halfedge.face(),
                                            barycoordsForHalfedgePoint(crossPoint.halfedge, crossPoint.tCross)));

    // Create a new face point in the begin face
    segmentPoints.push_front(
        SegmentEndpoint(crossPoint.halfedge.twin().face(),
                        barycoordsForHalfedgePoint(crossPoint.halfedge.twin(), 1.0 - crossPoint.tCross)));

  } else {

    // Remove the current ending face point
    segmentPoints.pop_front();

    // Remove the old crossing
    SegmentEndpoint crossPoint = segmentPoints.front();
    segmentPoints.pop_front();

    // Create a new face point in the face
    segmentPoints.push_front(
        SegmentEndpoint(crossPoint.halfedge.twin().face(),
                        barycoordsForHalfedgePoint(crossPoint.halfedge.twin(), 1.0 - crossPoint.tCross)));
  }
}

void MeshEmbeddedCurve::removeLastEndpoint() {

  if (segmentPoints.size() == 0) return;

  if (segmentPoints.size() == 1) {
    segmentPoints.pop_back();
    return;
  }

  // Break up a closed curve
  if (isClosed()) {

    // Remove the old crossing
    SegmentEndpoint crossPoint = segmentPoints.back();
    segmentPoints.pop_back();

    // Create a new face point in the end face
    segmentPoints.push_back(SegmentEndpoint(crossPoint.halfedge.face(),
                                            barycoordsForHalfedgePoint(crossPoint.halfedge, crossPoint.tCross)));

    // Create a new face point in the begin face
    segmentPoints.push_front(
        SegmentEndpoint(crossPoint.halfedge.twin().face(),
                        barycoordsForHalfedgePoint(crossPoint.halfedge.twin(), 1.0 - crossPoint.tCross)));

  } else {

    // Remove the current ending face point
    segmentPoints.pop_back();

    // Remove the old crossing
    SegmentEndpoint crossPoint = segmentPoints.back();
    segmentPoints.pop_back();

    // Create a new face point in the face
    segmentPoints.push_back(SegmentEndpoint(crossPoint.halfedge.face(),
                                            barycoordsForHalfedgePoint(crossPoint.halfedge, crossPoint.tCross)));
  }
}

void MeshEmbeddedCurve::rotateArbitraryStart() {

  if (!isClosed()) {
    throw std::runtime_error("Attempted to rotate non-closed curve");
  }

  SegmentEndpoint crossPoint = segmentPoints.front();
  segmentPoints.pop_front();
  segmentPoints.push_back(crossPoint);
}

void MeshEmbeddedCurve::closeCurve() {

  if (isClosed()) {
    throw std::runtime_error("Attempted to close curve which is already closed");
  }

  if (faceBefore(segmentPoints.front()) == faceAfter(segmentPoints.back())) {
    // Simply delete the two endpoints in the face!
    segmentPoints.pop_front();
    segmentPoints.pop_back();

  } else {
    // TODO handle general case by finding a shortest path or something?
    throw std::runtime_error("Tried to close curve for which endpoints do not lie in the same face");
  }
}

void MeshEmbeddedCurve::clearCurve() { segmentPoints.clear(); }


void MeshEmbeddedCurve::setFromZeroLevelset(VertexData<double>& implicitF) {

  // Note: If there are multiple disconnected zero sets, this finds any one of them abritrarily

  clearCurve();

  auto isForwardCrossingHalfedge = [&](HalfedgePtr he) {
    return (implicitF[he.vertex()] <= 0 && implicitF[he.twin().vertex()] > 0);
  };

  // = Find any halfedge crossing the from negative to positive
  HalfedgePtr startingHe;

  // Check boundary halfedges first
  for (HalfedgePtr he : mesh->realHalfedges()) {
    if (!he.twin().isReal()) {
      VertexPtr vTail = he.vertex();
      VertexPtr vTip = he.twin().vertex();

      if (isForwardCrossingHalfedge(he)) {
        startingHe = he;
        break;
      }
    }
  }
  // Check all halfedges now
  if (startingHe == HalfedgePtr()) {
    for (HalfedgePtr he : mesh->realHalfedges()) {
      VertexPtr vTail = he.vertex();
      VertexPtr vTip = he.twin().vertex();

      if (isForwardCrossingHalfedge(he)) {
        startingHe = he;
        break;
      }
    }
  }
  if (startingHe == HalfedgePtr()) {
    cout << "WARNING: Could not construct curve; implicit function has no zero level set" << endl;
    return;
  }

  // Add first point
  extendBack(startingHe.face(),
             barycoordsForHalfedgePoint(startingHe, scalarFunctionZeroPoint(implicitF[startingHe.vertex()],
                                                                            implicitF[startingHe.twin().vertex()])));

  // = Walk level set, building curve
  HalfedgePtr walkHe = startingHe;
  do {

    // Find next halfedge
    bool found = false;
    for (int i = 0; i < 2; i++) {
      walkHe = walkHe.next();

      if (!walkHe.edge().isBoundary() && isForwardCrossingHalfedge(walkHe.twin())) {
        found = true;
        break;
      }
    }

    // This should mean we hit a boundary
    if (!found) {
      break;
    }

    // Flip over the halfedge and add a point to the curve
    walkHe = walkHe.twin();

    extendBack(walkHe.face(),
               barycoordsForHalfedgePoint(
                   walkHe, scalarFunctionZeroPoint(implicitF[walkHe.vertex()], implicitF[walkHe.twin().vertex()])));


  } while (walkHe != startingHe);

  // Try to close the curve
  if (startingFace() == endingFace()) {
    closeCurve();
  }
}

FacePtr MeshEmbeddedCurve::startingFace() {
  if (isClosed() || segmentPoints.size() == 0) return FacePtr();
  return faceBefore(segmentPoints.front());
}


FacePtr MeshEmbeddedCurve::endingFace() {
  if (isClosed() || segmentPoints.size() == 0) return FacePtr();
  return faceAfter(segmentPoints.back());
}

std::vector<CurveSegment> MeshEmbeddedCurve::getCurveSegments() {

  std::vector<CurveSegment> segments;

  if (segmentPoints.size() == 0) return segments;

  for (size_t iP = 0; iP < segmentPoints.size() - 1; iP++) {

    // Construct the segment between this point and the next
    SegmentEndpoint& p1 = segmentPoints[iP];
    SegmentEndpoint& p2 = segmentPoints[iP + 1];

    CurveSegment newSeg;

    // Start point
    if (p1.isEdgeCrossing) {
      newSeg.face = p1.halfedge.twin().face();
      newSeg.startBaryCoord = barycoordsForHalfedgePoint(p1.halfedge.twin(), 1.0 - p1.tCross);
      newSeg.startPosition = positionOfSegmentEndpoint(p1);
      newSeg.startHe = p1.halfedge.twin();
    } else {
      newSeg.face = p1.face;
      newSeg.startBaryCoord = p1.faceCoords;
      newSeg.startPosition = positionOfSegmentEndpoint(p1);
      newSeg.startHe = HalfedgePtr();
    }

    // End point
    if (p2.isEdgeCrossing) {
      newSeg.endBaryCoord = barycoordsForHalfedgePoint(p2.halfedge, p2.tCross);
      newSeg.endPosition = positionOfSegmentEndpoint(p2);
      newSeg.endHe = p2.halfedge;
    } else {
      newSeg.endBaryCoord = p2.faceCoords;
      newSeg.endPosition = positionOfSegmentEndpoint(p2);
      newSeg.endHe = HalfedgePtr();
    }

    segments.push_back(newSeg);
  }

  // If the curve is closed, there is one more segment between the last and first entries
  if (isClosed()) {

    SegmentEndpoint& p1 = segmentPoints.back();
    SegmentEndpoint& p2 = segmentPoints.front();

    CurveSegment lastSeg;

    lastSeg.face = p1.halfedge.twin().face();
    lastSeg.startBaryCoord = barycoordsForHalfedgePoint(p1.halfedge.twin(), 1.0 - p1.tCross);
    lastSeg.startPosition = positionOfSegmentEndpoint(p1);
    lastSeg.startHe = p1.halfedge.twin();
    lastSeg.endBaryCoord = barycoordsForHalfedgePoint(p2.halfedge, p2.tCross);
    lastSeg.endPosition = positionOfSegmentEndpoint(p2);
    lastSeg.endHe = p2.halfedge;

    segments.push_back(lastSeg);
  }

  return segments;
}

bool MeshEmbeddedCurve::isClosed() {

  bool startIsClosed = segmentPoints.front().isEdgeCrossing;
  bool endIsClosed = segmentPoints.back().isEdgeCrossing;

  // Do a sanity check, because why not
  if (startIsClosed != endIsClosed) {
    throw std::runtime_error("Start and end of embedded curve disagree as to whether it is closed");
  }

  return startIsClosed;
}

void MeshEmbeddedCurve::validate() {

  // A segement of size 1 is silly. Allow zero, for newly constructed curves.
  if (segmentPoints.size() == 1) {
    throw std::runtime_error("MeshEmbeddedCurve segement point size == 1 doesn't make sense");
  }
  if (segmentPoints.size() == 0) {
    return;
  }

  // Either both or neither of the start and end should be closed
  bool startIsClosed = segmentPoints.front().isEdgeCrossing;
  bool endIsClosed = segmentPoints.back().isEdgeCrossing;
  if (startIsClosed != endIsClosed) {
    throw std::runtime_error("Start and end of embedded curve disagree as to whether it is closed");
  }

  // Check all interior segment points are not endpoints
  for (size_t iP = 1; iP < segmentPoints.size() - 1; iP++) {
    if (!segmentPoints[iP].isEdgeCrossing) {
      throw std::runtime_error("Interior points along embedded curve should not be endpoints");
    }
  }

  // Check that the path is well-connected
  for (size_t iP = 0; iP < segmentPoints.size() - 1; iP++) {
    SegmentEndpoint& p1 = segmentPoints[iP];
    SegmentEndpoint& p2 = segmentPoints[iP + 1];

    if (faceAfter(p1) != faceBefore(p2)) {
      throw std::runtime_error("Embedded curve path segment points do not describe a path through faces");
    }
  }
  if (isClosed()) {
    SegmentEndpoint& p1 = segmentPoints.back();
    SegmentEndpoint& p2 = segmentPoints.front();
    if (faceAfter(p1) != faceBefore(p2)) {
      throw std::runtime_error("Embedded curve path segment points do not describe a path through faces");
    }
  }
}

double MeshEmbeddedCurve::computeLength() {
  double l = 0;
  for (CurveSegment& c : getCurveSegments()) {
    l += c.length();
  }
  return l;
}


size_t MeshEmbeddedCurve::nSegments() {
  if(isClosed()) {
    return segmentPoints.size();
  } else {
    return segmentPoints.size()-1;
  }
}

double CurveSegment::length() { return norm(startPosition - endPosition); }


MeshEmbeddedCurve MeshEmbeddedCurve::copy(HalfedgeMeshDataTransfer& transfer, Geometry<Euclidean>* otherGeom) {

  MeshEmbeddedCurve newCurve(otherGeom);

  // Copy each segment
  for (SegmentEndpoint& e : segmentPoints) {
    if (e.isEdgeCrossing) {
      newCurve.segmentPoints.push_back(SegmentEndpoint(transfer.heMapBack[e.halfedge], e.tCross));
    } else {
      newCurve.segmentPoints.push_back(SegmentEndpoint(transfer.fMapBack[e.face], e.faceCoords));
    }
  }

  return newCurve;
}

MeshEmbeddedCurve MeshEmbeddedCurve::copyBack(HalfedgeMeshDataTransfer& transfer, Geometry<Euclidean>* otherGeom) {

  MeshEmbeddedCurve newCurve(otherGeom);

  // Copy each segment
  for (SegmentEndpoint& e : segmentPoints) {
    if (e.isEdgeCrossing) {
      newCurve.segmentPoints.push_back(SegmentEndpoint(transfer.heMap[e.halfedge], e.tCross));
    } else {
      newCurve.segmentPoints.push_back(SegmentEndpoint(transfer.fMap[e.face], e.faceCoords));
    }
  }

  return newCurve;
}
