#include "NormalCoordinateIntrinsicTriangulation.h"

NormalCoordinateIntrinsicTriangulation::NormalCoordinateIntrinsicTriangulation(
    ManifoldSurfaceMesh& mesh_, IntrinsicGeometryInterface& inputGeom_)
    // Note: this initializer list does something slightly wacky: it creates the
    // new mesh on the heap, then loses track of pointer while setting the
    // BaseGeometryInterface::mesh reference to it. Later, it picks the pointer
    // back up from the reference and wraps it in the intrinsicMesh
    // unique_ptr<>. I believe that this is all valid, but its probably a sign
    // of bad design.
    : IntrinsicGeometryInterface(*mesh_.copy().release()), inputMesh(mesh_),
      inputGeom(inputGeom_),
      intrinsicMesh(dynamic_cast<ManifoldSurfaceMesh*>(&mesh)),
      normalCoordinates_ptr(
          NormalCoordinates::constructFromEdges(*intrinsicMesh)),
      normalCoordinates(*normalCoordinates_ptr) {

    // Make sure the input mesh is triangular
    if (!mesh.isTriangular()) {
        throw std::runtime_error(
            "normal coordinate triangulation requires triangle mesh as input");
    }

    inputGeom.requireEdgeLengths();
    // normalCoordinates =
    // NormalCoordinates::constructFromEdges(*intrinsicMesh);
    intrinsicEdgeLengths = inputGeom.edgeLengths.reinterpretTo(mesh);
    inputGeom.unrequireEdgeLengths();
    // TODO: if I don't unrequire this, then inserting vertices is really slow
    // for some reason. But I do want to keep this buffer around - I just know
    // it should never update

    mollifyIntrinsic(*intrinsicMesh, intrinsicEdgeLengths, 1e-5);

    edgeLengths = EdgeData<double>(*intrinsicMesh);
    for (Edge e : intrinsicMesh->edges())
        edgeLengths[e] = intrinsicEdgeLengths[e];

    intrinsicToInput = VertexData<SurfacePoint>(*intrinsicMesh);
    inputToIntrinsic = VertexData<SurfacePoint>(inputMesh);
    for (size_t iV = 0; iV < intrinsicMesh->nVertices(); iV++) {
        Vertex vIntrinsic            = intrinsicMesh->vertex(iV);
        Vertex vInput                = inputMesh.vertex(iV);
        intrinsicToInput[vIntrinsic] = vInput;
        inputToIntrinsic[vInput]     = vIntrinsic;
    }

    // We will manually keep these buffers up to date while mutating the mesh
    // This is a crucial performance optimization
    // (Also we'll manually update edge lengths)
    requireVertexAngleSums();
    requireCornerAngles();
    requireHalfedgeVectorsInVertex();
    requireHalfedgeVectorsInFace();

    // Boundary edges are always fixed
    isFixed = EdgeData<bool>(*intrinsicMesh, false);
    for (BoundaryLoop b : intrinsicMesh->boundaryLoops()) {
        for (Edge e : b.adjacentEdges()) {
            isFixed[e] = true;
            // fixedEdges.push_back(e);
        }
    }
}

// ======================================================
//                 Queries & Accesses
// ======================================================

// TODO
CommonSubdivision&
NormalCoordinateIntrinsicTriangulation::traceEdges(bool geodesic) {

    if (!commonSubdivision) {
        commonSubdivision =
            std::make_unique<CommonSubdivision>(inputMesh, *intrinsicMesh);
        traceEdges(*commonSubdivision, geodesic);
    }
    return *commonSubdivision;
}


// Return by reference to make tests easier to write
void NormalCoordinateIntrinsicTriangulation::traceEdges(CommonSubdivision& cs,
                                                        bool geodesic) {
    // Mesh A is inputMesh, Mesh B is intrinsicMesh
    GC_SAFETY_ASSERT(&cs.meshA == &inputMesh,
                     "CommonSubdivision.meshA must be inputMesh");
    GC_SAFETY_ASSERT(&cs.meshB == &(*intrinsicMesh),
                     "CommonSubdivision.meshB must be intrinsicMesh");

    verbose_assert(inputMesh.nVertices() <= intrinsicMesh->nVertices(),
                   "tracing simplified edges not implemented yet");

    intrinsicMesh->compress();


    // === Debugging info
    if (false) {
        for (Face f : intrinsicMesh->faces()) {
            bool interesting = false;
            for (Vertex v : f.adjacentVertices()) {
                if (intrinsicToInput[v].type != SurfacePointType::Vertex) {
                    interesting = true;
                    break;
                }
            }
            if (interesting) {
                cout << endl << endl;
                cout << "INTERESTING FACE " << f << myendl;
                for (Vertex v : f.adjacentVertices()) {
                    cout << "\t\t Vertex " << intrinsicToInput[v] << myendl;
                }
            }
        }
    }


    // TODO: Don't require meshB's vertices to be a superset of meshA's

    // Construct CommonSubdivisionPoints corresponding to shared vertices
    VertexData<CommonSubdivisionPoint*> aVtx(inputMesh);
    VertexData<CommonSubdivisionPoint*> bVtx(*intrinsicMesh);

    for (size_t iV = 0; iV < intrinsicMesh->nVertices(); iV++) {
        Vertex vB       = intrinsicMesh->vertex(iV);
        SurfacePoint pA = intrinsicToInput[vB];

        switch (pA.type) {
        case SurfacePointType::Vertex: {
            Vertex vA = pA.vertex;

            cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                CSIntersectionType::VERTEX_VERTEX, pA, SurfacePoint(vB), true});

            aVtx[vA] = &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];
            break;
        }
        case SurfacePointType::Edge:
            cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                CSIntersectionType::EDGE_VERTEX, pA, SurfacePoint(vB), true});
            break;
        case SurfacePointType::Face:
            cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                CSIntersectionType::FACE_VERTEX, pA, SurfacePoint(vB), true});
            break;
        }

        bVtx[vB] = &cs.subdivisionPoints.back();
    }

    // Allocate space for points along B's edges, and also fill in the endpoints
    // of B's edges
    // Note that for shared edges, we will actually need 3 entries instead of 2
    // (due to a HACK), but we'll fix that later
    for (Edge eB : intrinsicMesh->edges()) {
        int n = positivePart(normalCoordinates[eB]);
        cs.pointsAlongB[eB].resize(n + 2);
        cs.pointsAlongB[eB][0]     = bVtx[src(eB)];
        cs.pointsAlongB[eB][n + 1] = bVtx[dst(eB)];
    }

    // Trace the edges of mesh A (inputMesh) over mesh B (intrinsicMesh)
    for (Edge eA : inputMesh.edges()) {
        CompoundCurve compoundPath = traceInputEdge(eA);

        for (size_t iC = 0; iC < compoundPath.components.size(); iC++) {

            const Curve& curve = compoundPath.components[iC];
            bool first         = iC == 0;

            auto& path = curve.crossings;

            if (path.size() == 1 && std::get<0>(path[0]) < 0) {

                // Shared edge
                Halfedge heB             = std::get<1>(path[0]);
                Edge eB                  = heB.edge();
                bool positiveOrientation = heB == eB.halfedge();

                verbose_assert(normalCoordinates[eB] < 0,
                               "eB should have negative n, but has n = " +
                                   std::to_string(normalCoordinates[eB]));

                // make pointsAlongB[eB] one longer
                cs.pointsAlongB[eB].push_back(
                    cs.pointsAlongB[eB][cs.pointsAlongB[eB].size() - 1]);

                // Construct intersection point, flagging it with barycentric
                // coordinate 0.5 to indicate that the edges are parallel
                cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                    CSIntersectionType::EDGE_PARALLEL, SurfacePoint(eA, 0.5),
                    SurfacePoint(heB.edge(), 0.5), positiveOrientation});

                CommonSubdivisionPoint* crPt =
                    &cs.subdivisionPoints[cs.subdivisionPoints.size() - 1];

                if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);
                cs.pointsAlongA[eA].push_back(crPt);
                cs.pointsAlongA[eA].push_back(bVtx[dst(eB)]);

                cs.pointsAlongB[eB][1] = crPt;
            } else {
                // Indirect path

                if (geodesic) { // Lay out geodesic
                    // cout << path << myendl;

                    std::vector<std::pair<SurfacePoint, double>> geodesicPath =
                        generateFullSingleGeodesicGeometry(*intrinsicMesh,
                                                           *this, curve);

                    if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);

                    Halfedge hB;
                    for (size_t iC = 0; iC < path.size(); ++iC) {
                        hB                       = std::get<1>(path[iC]);
                        Edge eB                  = hB.edge();
                        bool positiveOrientation = hB == eB.halfedge();

                        // geodesicPath stores the start and end point, which
                        // path doesn't do, so we offset by 1 here
                        SurfacePoint ptB = std::get<0>(geodesicPath[iC + 1]);
                        double tA        = std::get<1>(geodesicPath[iC + 1]);

                        int iB = std::get<0>(path[iC]);
                        if (!positiveOrientation)
                            iB = positivePart(normalCoordinates[eB]) - iB - 1;

                        cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                            CSIntersectionType::EDGE_TRANSVERSE,
                            SurfacePoint(eA, tA), ptB, positiveOrientation});

                        CommonSubdivisionPoint* crPt =
                            &cs.subdivisionPoints[cs.subdivisionPoints.size() -
                                                  1];

                        cs.pointsAlongA[eA].push_back(crPt);
                        cs.pointsAlongB[eB][iB + 1] = crPt;
                    }

                    Vertex bDst = hB.twin().next().tipVertex();
                    cs.pointsAlongA[eA].push_back(bVtx[bDst]);
                } else { // Space points evenly

                    // cout << "PATH" << endl;

                    if (first) cs.pointsAlongA[eA].push_back(aVtx[src(eA)]);
                    Halfedge hB;
                    for (const std::pair<int, Halfedge>& pt : path) {
                        hB                       = std::get<1>(pt);
                        Edge eB                  = hB.edge();
                        bool positiveOrientation = hB == eB.halfedge();

                        int iB = std::get<0>(pt);
                        if (!positiveOrientation)
                            iB = positivePart(normalCoordinates[eB]) - iB - 1;

                        double tB =
                            (iB + 1) /
                            (double)(positivePart(normalCoordinates[eB]) + 1);

                        cs.subdivisionPoints.push_back(CommonSubdivisionPoint{
                            CSIntersectionType::EDGE_TRANSVERSE,
                            SurfacePoint(eA, 0.5), SurfacePoint(hB.edge(), tB),
                            positiveOrientation});

                        CommonSubdivisionPoint* crPt =
                            &cs.subdivisionPoints[cs.subdivisionPoints.size() -
                                                  1];

                        cs.pointsAlongA[eA].push_back(crPt);
                        cs.pointsAlongB[eB][iB + 1] = crPt;
                    }
                    Vertex bDst = hB.twin().next().tipVertex();
                    cs.pointsAlongA[eA].push_back(bVtx[bDst]);

                    // Space out points on mesh A
                    for (size_t iA = 1; iA + 1 < cs.pointsAlongA[eA].size();
                         iA++) {
                        double tA =
                            iA / (double)(cs.pointsAlongA[eA].size() - 1);
                        cs.pointsAlongA[eA][iA]->posA.tEdge = tA;
                    }
                }
            }
        }
    }

    // Check that we've accounted for all promised crossings along B's edges
    for (Edge eB : intrinsicMesh->edges()) {
        for (const auto& p : cs.pointsAlongB[eB]) {
            if (!p) {
                cout << "Missed a crossing on edge " << eB
                     << " which is supposed to have " << normalCoordinates[eB]
                     << myendl;
                cout << "Endpoints "
                     << intrinsicToInput[eB.halfedge().tailVertex()] << myendl
                     << "\t->" << intrinsicToInput[eB.halfedge().tipVertex()]
                     << myendl;
                verbose_assert(p, "oops, missed a crossing");
            }
        }
        for (size_t iC = 1; iC + 1 < cs.pointsAlongB[eB].size(); ++iC) {
            if (cs.pointsAlongB[eB][iC]->intersectionType ==
                CSIntersectionType::VERTEX_VERTEX) {
                cout << "Error at crossing " << iC << " of "
                     << cs.pointsAlongB[eB].size() << " on edge " << eB
                     << myendl;
                for (const auto& pt : cs.pointsAlongB[eB])
                    cout << "\t" << *pt << myendl;
                throw_verbose_runtime_error("encountered vertex intersection "
                                            "in the middle of an edge ?!");
            }
        }

        // SurfacePoint src = intrinsicToInput[eB.halfedge().tailVertex()];
        // SurfacePoint dst = intrinsicToInput[eB.halfedge().tipVertex()];
        // if (src.type != SurfacePointType::Vertex ||
        //     dst.type != SurfacePointType::Vertex) {
        //     cout << endl << "Interesting edge " << eB << myendl;
        //     WATCH2(src, dst);
        //     for (const auto& p : cs.pointsAlongB[eB]) {
        //         cout << "\t" << *p << myendl;
        //     }
        // }
    }
}

std::vector<std::vector<SurfacePoint>>
NormalCoordinateIntrinsicTriangulation::traceEdgeSet() {

    CommonSubdivision& subd = traceEdges();

    std::vector<std::vector<SurfacePoint>> out;
    for (Edge e : intrinsicMesh->edges()) {
        std::vector<SurfacePoint> thisEdge;
        for (CommonSubdivisionPoint* p : subd.pointsAlongB[e]) {
            thisEdge.push_back(p->posA);
        }
        out.push_back(thisEdge);
    }


    return out;
}

// ======================================================
//                High-Level Mutators
// ======================================================

size_t NormalCoordinateIntrinsicTriangulation::flipToDelaunay(double tol) {
    std::deque<Edge> edgesToCheck;
    EdgeData<bool> inQueue(mesh, true);
    for (Edge e : mesh.edges()) {
        edgesToCheck.push_back(e);
    }

    size_t nFlips = 0;
    while (!edgesToCheck.empty()) {

        // Get the top element from the queue of possibily non-Delaunay
        // edges
        Edge e = edgesToCheck.front();
        edgesToCheck.pop_front();
        inQueue[e] = false;

        bool wasFlipped = flipEdgeIfNotDelaunay(e, tol);
        if (e.getIndex() == 4849 && false) {
            cout << " Considering edge " << e << ": wasFlipped: " << wasFlipped
                 << "\t isDelaunay: " << isDelaunay(e, tol) << myendl;
            std::array<Halfedge, 6> diamond{e.halfedge(),
                                            e.halfedge().next(),
                                            e.halfedge().next().next(),
                                            e.halfedge().twin(),
                                            e.halfedge().twin().next(),
                                            e.halfedge().twin().next().next()};
            cout << " Diamond edge lengths: " << myendl;
            for (Halfedge he : diamond) {
                cout << "\t\t " << intrinsicEdgeLengths[he.edge()] << myendl;
            }
            cout << " Diamond corner angles: " << myendl;
            for (Halfedge he : diamond) {
                cout << "\t\t " << (cornerAngles[he.corner()] / M_PI * 180.)
                     << " degrees" << myendl;
            }
        }

        if (!wasFlipped) continue;

        // Handle the aftermath of a flip
        nFlips++;

        // Add neighbors to queue, as they may need flipping now
        Halfedge he                  = e.halfedge();
        Halfedge heN                 = he.next();
        Halfedge heT                 = he.twin();
        Halfedge heTN                = heT.next();
        std::vector<Edge> neighEdges = {heN.edge(), heN.next().edge(),
                                        heTN.edge(), heTN.next().edge()};
        for (Edge nE : neighEdges) {
            if (!inQueue[nE]) {
                edgesToCheck.push_back(nE);
                inQueue[nE] = true;
            }
        }
    }

    for (Edge e : mesh.edges()) {
        if (!isDelaunay(e, tol)) {
            WATCH2(getEdgeCotanWeight(e), e.isBoundary());
        }
        verbose_assert(isDelaunay(e, tol),
                       "edge " + std::to_string(e) + " not delaunay?!");
    }

    // refreshQuantities(); // update buffers manually instead
    return nFlips;
}

// ======================================================
//                Low-Level Mutators
// ======================================================

// If the edge is not Delaunay, flip it. Returns true if flipped.
// TODO
bool NormalCoordinateIntrinsicTriangulation::flipEdgeIfNotDelaunay(
    Edge e, double possibleEPS) {
    // Can't flip
    if (isFixed[e]) return false;

    if (isDelaunay(e, possibleEPS)) return false;

    // Get geometric data
    Halfedge he                            = e.halfedge();
    std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

    /*
    // Test if geometryically flippable flippable (both signed areas of new
    // triangles are positive)
    double A1      = cross(layoutPositions[1] - layoutPositions[0],
                      layoutPositions[3] - layoutPositions[0]);
    double A2      = cross(layoutPositions[3] - layoutPositions[2],
                      layoutPositions[1] - layoutPositions[2]);
    double areaEPS = possibleEPS * (A1 + A2);
    if (A1 < areaEPS || A2 < areaEPS) {
        return false;
    }
    */


    // Compute the new edge length
    double newLength = (layoutPositions[1] - layoutPositions[3]).norm();


    // Combinatorial flip
    auto nUpdate = normalCoordinates.computeFlippedData(e);
    bool flipped = intrinsicMesh->flip(e, false);

    // Might not have been flippable for connectivity reasons
    if (!flipped) {
        return false;
    }

    // If we're going to create a non-finite edge length, abort the flip
    // (only happens if you're in a bad numerical place)
    if (!std::isfinite(newLength)) {
        intrinsicMesh->flip(e, false);
        return false;
    }

    // Assign the new edge lengths
    // TODO project to satisfy triangle inequality?
    intrinsicEdgeLengths[e] = newLength;
    edgeLengths[e]          = newLength;
    normalCoordinates.applyFlippedData(e, nUpdate);

    // === Update various quantities

    // depends on intrinsicEdgeLengths
    updateFaceArea(e.halfedge().face());
    updateFaceArea(e.halfedge().twin().face());

    // depends on intrinsicEdgeLengths, faceAreas
    updateHalfedgeVectorsInFace(e.halfedge().face());
    updateHalfedgeVectorsInFace(e.halfedge().twin().face());

    std::array<Corner, 6> incidentCorners{
        e.halfedge().corner(),
        e.halfedge().next().corner(),
        e.halfedge().next().next().corner(),
        e.halfedge().twin().corner(),
        e.halfedge().twin().next().corner(),
        e.halfedge().twin().next().next().corner()};
    for (Corner c : incidentCorners) {
        // depends on intrinsicEdgeLengths
        updateCornerAngle(c);
    }

    std::array<Vertex, 4> incidentVertices{
        e.halfedge().vertex(), e.halfedge().next().vertex(),
        e.halfedge().next().next().vertex(),
        e.halfedge().twin().next().next().vertex()};
    for (Vertex v : incidentVertices) {
        // depends on cornerAngles
        updateVertexAngleSum(v);

        // depends on cornerAngles, vertexAngleSums
        updateHalfedgeVectorsInVertex(v);
    }

    // Do callbacks
    invokeEdgeFlipCallbacks(e);

    return true;
}

// If the edge can be flipped, flip it (must be combinatorially flippable
// and inside a convex quad). Returns true if flipped.
// TODO
bool NormalCoordinateIntrinsicTriangulation::flipEdgeIfPossible(
    Edge e, double possibleEPS, bool verbose) {
    // Can't flip
    if (isFixed[e]) return false;

    // Get geometric data
    Halfedge he                            = e.halfedge();
    std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

    // Test if geometryically flippable flippable (both signed areas of new
    // triangles are positive)
    double A1      = cross(layoutPositions[1] - layoutPositions[0],
                      layoutPositions[3] - layoutPositions[0]);
    double A2      = cross(layoutPositions[3] - layoutPositions[2],
                      layoutPositions[1] - layoutPositions[2]);
    double areaEPS = possibleEPS * (A1 + A2);

    if (A1 < areaEPS || A2 < areaEPS) {
        // std::cout << "  degen eps = "
        //           << std::min((A1 - areaEPS), (A2 - areaEPS)) << std::endl;
        return false;
    }


    // Combinatorial flip
    auto nUpdate = normalCoordinates.computeFlippedData(e);
    bool flipped = intrinsicMesh->flip(e, false);

    // Might not have been flippable for connectivity reasons
    if (!flipped) {
        return false;
    }

    // Compute the new edge length
    double newLength = (layoutPositions[1] - layoutPositions[3]).norm();

    // If we're going to create a non-finite edge length, abort the flip
    // (only happens if you're in a bad numerical place)
    if (!std::isfinite(newLength)) {
        intrinsicMesh->flip(e, false);

        return false;
    }

    // Assign the new edge lengths
    // TODO project to satisfy triangle inequality?
    intrinsicEdgeLengths[e] = newLength;
    edgeLengths[e]          = newLength;
    normalCoordinates.applyFlippedData(e, nUpdate);

    // === Update various quantities

    // depends on intrinsicEdgeLengths
    updateFaceArea(e.halfedge().face());
    updateFaceArea(e.halfedge().twin().face());

    // depends on intrinsicEdgeLengths, faceAreas
    updateHalfedgeVectorsInFace(e.halfedge().face());
    updateHalfedgeVectorsInFace(e.halfedge().twin().face());

    std::array<Corner, 6> incidentCorners{
        e.halfedge().corner(),
        e.halfedge().next().corner(),
        e.halfedge().next().next().corner(),
        e.halfedge().twin().corner(),
        e.halfedge().twin().next().corner(),
        e.halfedge().twin().next().next().corner()};
    for (Corner c : incidentCorners) {
        // depends on intrinsicEdgeLengths
        updateCornerAngle(c);
    }

    std::array<Vertex, 4> incidentVertices{
        e.halfedge().vertex(), e.halfedge().next().vertex(),
        e.halfedge().next().next().vertex(),
        e.halfedge().twin().next().next().vertex()};
    for (Vertex v : incidentVertices) {
        // depends on cornerAngles
        updateVertexAngleSum(v);

        // depends on cornerAngles, vertexAngleSums
        updateHalfedgeVectorsInVertex(v);
    }


    // Do callbacks
    invokeEdgeFlipCallbacks(e);

    return true;
}

double NormalCoordinateIntrinsicTriangulation::checkFlip(Edge e) {
    // Can't flip
    if (isFixed[e]) return std::numeric_limits<double>::infinity();

    // Check topologically flippable
    {
        Halfedge ha1 = e.halfedge();
        Halfedge ha2 = ha1.next();
        Halfedge ha3 = ha2.next();
        Halfedge hb1 = ha1.sibling();
        Halfedge hb2 = hb1.next();
        Halfedge hb3 = hb2.next();

        // incident on degree 1 vertex
        if (ha2 == hb1 || hb2 == ha1) {
            return std::numeric_limits<double>::infinity();
        }
    }

    // Get geometric data
    Halfedge he                            = e.halfedge();
    std::array<Vector2, 4> layoutPositions = layoutDiamond(he);

    // Test if geometryically flippable flippable (both signed areas of new
    // triangles are positive)
    double A1      = cross(layoutPositions[1] - layoutPositions[0],
                      layoutPositions[3] - layoutPositions[0]);
    double A2      = cross(layoutPositions[3] - layoutPositions[2],
                      layoutPositions[1] - layoutPositions[2]);
    double areaSum = (A1 + A2);

    return std::min(A1 / areaSum, A2 / areaSum);
}

Vertex NormalCoordinateIntrinsicTriangulation::insertCircumcenterOrSplitSegment(
    Face f, bool verbose) {
    // === Circumcenter in barycentric coordinates

    Halfedge he0            = f.halfedge();
    double a                = intrinsicEdgeLengths[he0.next().edge()];
    double b                = intrinsicEdgeLengths[he0.next().next().edge()];
    double c                = intrinsicEdgeLengths[he0.edge()];
    double a2               = a * a;
    double b2               = b * b;
    double c2               = c * c;
    Vector3 circumcenterLoc = {a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2),
                               c2 * (a2 + b2 - c2)};
    circumcenterLoc         = normalizeBarycentric(circumcenterLoc);

    // Trace from the barycenter (have to trace from somewhere)
    Vector3 barycenter        = Vector3::constant(1. / 3.);
    Vector3 vecToCircumcenter = circumcenterLoc - barycenter;

    // === Trace the ray to find the location of the new point on the intrinsic
    // meshes

    // Data we need from the intrinsic trace
    TraceOptions options;
    // if (markedEdges.size() > 0) { // TODO: what is this?
    //     options.barrierEdges = &markedEdges;
    // }
    TraceGeodesicResult intrinsicTraceResult =
        traceGeodesic(*this, f, barycenter, vecToCircumcenter, options);
    SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;
    if (newPositionOnIntrinsic.type == SurfacePointType::Vertex &&
        newPositionOnIntrinsic.vertex == Vertex()) {
        // tracing failed
        return Vertex();
    }

    // If the circumcenter is blocked by an edge, insert the midpoint of that
    // edge instead (which happens to be just what is needed for Chew's 2nd
    // algo).
    // TODO: also do this if we hit locked edges
    // TODO: delete inserted vertices inside the diametral circle?
    if (newPositionOnIntrinsic.type == SurfacePointType::Edge) {
        newPositionOnIntrinsic.tEdge = 0.5;
    }

    if (verbose) {
        cout << "insertion point: " << newPositionOnIntrinsic << myendl;
    }

    // === Phase 3: Add the new vertex
    Vertex newV = insertVertex(newPositionOnIntrinsic, verbose);
    return newV;
}

// Assumes intrinsicEdgeLengths is up to date
void NormalCoordinateIntrinsicTriangulation::updateCornerAngle(Corner c) {
    cornerAngles[c] = getCornerAngle(c);
}

// Assumes cornerAngles exists and is up to date
void NormalCoordinateIntrinsicTriangulation::updateVertexAngleSum(Vertex v) {
    vertexAngleSums[v] = 0;
    for (Corner c : v.adjacentCorners()) vertexAngleSums[v] += cornerAngles[c];
}

// Assumes that intrinsicEdgeLengths is up to date
void NormalCoordinateIntrinsicTriangulation::updateFaceArea(Face f) {
    faceAreas[f] = getFaceArea(f);
}

// Assumes cornerAngles, vertexAngleSums exist and are up to date
void NormalCoordinateIntrinsicTriangulation::updateHalfedgeVectorsInVertex(
    Vertex v) {
    // stolen from intrinsic_geometry_interface.cpp

    auto cornerScaledAngle = [&](Corner c) -> double {
        if (c.vertex().isBoundary()) {
            double s = PI / vertexAngleSums[c.vertex()];
            return s * cornerAngles[c];
        } else {
            double s = 2.0 * PI / vertexAngleSums[c.vertex()];
            return s * cornerAngles[c];
        }
    };

    double coordSum = 0.0;

    // Custom loop to orbit CCW
    Halfedge firstHe = v.halfedge();
    Halfedge currHe  = firstHe;
    do {
        halfedgeVectorsInVertex[currHe] =
            Vector2::fromAngle(coordSum) * edgeLengths[currHe.edge()];
        coordSum += cornerScaledAngle(currHe.corner());
        if (!currHe.isInterior()) break;
        currHe = currHe.next().next().twin();
    } while (currHe != firstHe);
}

// Assumes intrinsicEdgeLengths, faceAreas exist and are up to date
void NormalCoordinateIntrinsicTriangulation::updateHalfedgeVectorsInFace(
    Face f) {
    // stolen from intrinsic_geometry_interface.cpp

    // Gather some values
    Halfedge heAB = f.halfedge();
    Halfedge heBC = heAB.next();
    Halfedge heCA = heBC.next();
    GC_SAFETY_ASSERT(heCA.next() == heAB, "faces must be triangular");

    double lAB = intrinsicEdgeLengths[heAB.edge()];
    double lBC = intrinsicEdgeLengths[heBC.edge()];
    double lCA = intrinsicEdgeLengths[heCA.edge()];

    // Assign positions to all three vertices
    // Vector2 pA{0., 0.}; // used implicitly
    Vector2 pB{lAB, 0.};
    // pC is the hard one:

    double tArea = faceAreas[f];

    // Compute width and height of right triangle formed via altitude from C
    double h = 2. * tArea / lAB;
    double w = std::sqrt(std::max(0., lCA * lCA - h * h));

    // Take the closer of the positive and negative solutions
    if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;

    // Project some vectors to get the actual position
    Vector2 pC{w, h};

    // Now, all halfedge vectors are just coordinates
    halfedgeVectorsInFace[heAB] = pB;
    halfedgeVectorsInFace[heBC] = pC - pB;
    halfedgeVectorsInFace[heCA] = -pC;
}

// ======================================================
//                Low-Level Queries
// ======================================================
double NormalCoordinateIntrinsicTriangulation::getCornerAngle(Corner c) const {
    double angle = std::acos(getCornerAngleCosine(c));

    return angle;
}
double
NormalCoordinateIntrinsicTriangulation::getCornerAngleCosine(Corner c) const {
    Halfedge heA   = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB   = heOpp.next();

    GC_SAFETY_ASSERT(heB.next() == heA, "faces mush be triangular");

    double lOpp = intrinsicEdgeLengths[heOpp.edge()];
    double lA   = intrinsicEdgeLengths[heA.edge()];
    double lB   = intrinsicEdgeLengths[heB.edge()];

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
    q        = clamp(q, -1.0, 1.0);

    return q;
}

double NormalCoordinateIntrinsicTriangulation::getFaceArea(Face f) const {
    Halfedge he = f.halfedge();
    double a    = intrinsicEdgeLengths[he.edge()];
    he          = he.next();
    double b    = intrinsicEdgeLengths[he.edge()];
    he          = he.next();
    double c    = intrinsicEdgeLengths[he.edge()];
    return triangleArea(a, b, c);
}

double NormalCoordinateIntrinsicTriangulation::getHalfedgeCotanWeight(
    Halfedge heI) const {
    if (heI.isInterior()) {
        Halfedge he = heI;
        double l_ij = intrinsicEdgeLengths[he.edge()];
        he          = he.next();
        double l_jk = intrinsicEdgeLengths[he.edge()];
        he          = he.next();
        double l_ki = intrinsicEdgeLengths[he.edge()];
        he          = he.next();
        GC_SAFETY_ASSERT(he == heI, "faces mush be triangular");
        return halfedgeCotanWeight(l_ij, l_jk, l_ki);
    } else {
        return 0.;
    }
}

double
NormalCoordinateIntrinsicTriangulation::getEdgeCotanWeight(Edge e) const {
    return getHalfedgeCotanWeight(e.halfedge()) +
           getHalfedgeCotanWeight(e.halfedge().twin());
}

double NormalCoordinateIntrinsicTriangulation::getCircumradius(Face f) const {
    double a = intrinsicEdgeLengths[f.halfedge().edge()];
    double b = intrinsicEdgeLengths[f.halfedge().next().edge()];
    double c = intrinsicEdgeLengths[f.halfedge().next().next().edge()];

    double area = triangleArea(a, b, c);

    return a * b * c / (4.0 * area);
}

double NormalCoordinateIntrinsicTriangulation::getShortestEdge(Face f) const {
    double a = intrinsicEdgeLengths[f.halfedge().edge()];
    double b = intrinsicEdgeLengths[f.halfedge().next().edge()];
    double c = intrinsicEdgeLengths[f.halfedge().next().next().edge()];

    return fmin(a, fmin(b, c));
}

double NormalCoordinateIntrinsicTriangulation::getMinAngle() const {
    double minCornerAngle = 10;
    for (Corner c : intrinsicMesh->corners()) {
        minCornerAngle = fmin(minCornerAngle, cornerAngles[c]);
    }
    return minCornerAngle;
}

double
NormalCoordinateIntrinsicTriangulation::getMaxCircumcircleToEdgeLengthRatio()
    const {
    double maxRatio = -1;
    for (Face f : intrinsicMesh->faces()) {
        double cr = getCircumradius(f);
        double el = getShortestEdge(f);
        maxRatio  = fmax(maxRatio, cr / el);
    }
    return maxRatio;
}

double NormalCoordinateIntrinsicTriangulation::getMinInteriorAngleSum() const {
    double minSum = std::numeric_limits<double>::infinity();
    for (Vertex v : intrinsicMesh->vertices()) {
        if (!v.isBoundary()) minSum = fmin(vertexAngleSums[v], minSum);
    }
    return minSum;
}

double NormalCoordinateIntrinsicTriangulation::getMinBoundaryAngleSum() const {
    double minSum = std::numeric_limits<double>::infinity();
    for (Vertex v : intrinsicMesh->vertices()) {
        if (v.isBoundary()) minSum = fmin(vertexAngleSums[v], minSum);
    }
    return minSum;
}

// Get min angle on faces whose vertices all have angle sum at least
// minAngleSum
double NormalCoordinateIntrinsicTriangulation::getMinAngleOnValidFaces(
    double minAngleSum) const {

    double minCornerAngle = 10;
    for (Face f : intrinsicMesh->faces()) {
        if (!faceHasLargeAngleSums(f, minAngleSum)) continue;
        if (!parentFaceHasLargeAngleSums(f, minAngleSum)) continue;
        for (Corner c : f.adjacentCorners()) {
            minCornerAngle = fmin(minCornerAngle, cornerAngles[c]);
        }
    }
    return minCornerAngle;
}

// Get max cimrcumcircle-to-edge-length ratio on faces whose vertices all
// have angle sum at least minAngleSum
double NormalCoordinateIntrinsicTriangulation::getMaxRatioOnValidFaces(
    double minAngleSum) const {

    double maxRatio = -1;
    for (Face f : intrinsicMesh->faces()) {
        if (!faceHasLargeAngleSums(f, minAngleSum)) continue;
        if (!parentFaceHasLargeAngleSums(f, minAngleSum)) continue;

        double cr = getCircumradius(f);
        double el = getShortestEdge(f);
        maxRatio  = fmax(maxRatio, cr / el);
    }
    return maxRatio;
}

// Check if all face vertices have angle sum at least minAngleSum
bool NormalCoordinateIntrinsicTriangulation::faceHasLargeAngleSums(
    Face f, double minAngleSum) const {
    for (Vertex v : f.adjacentVertices()) {
        if (vertexAngleSums[v] < M_PI / 3) {
            return false;
        }
    }
    return true;
}

// Check if all face vertices have angle sum at least minAngleSum
bool NormalCoordinateIntrinsicTriangulation::parentFaceHasLargeAngleSums(
    Face f, double minAngleSum) const {
    Face fInput = getParentFace(f);
    if (fInput == Face()) return true;

    inputGeom.requireVertexAngleSums();
    for (Vertex v : fInput.adjacentVertices()) {
        if (inputGeom.vertexAngleSums[v] < M_PI / 3) {
            return false;
        }
    }
    inputGeom.unrequireVertexAngleSums();
    return true;
}

bool NormalCoordinateIntrinsicTriangulation::isOnFixedEdge(Vertex v) const {
    for (Edge e : v.adjacentEdges()) {
        if (isFixed[e]) return true;
    }
    return false;
}

// Takes in a halfedge of the intrinsic mesh whose edge's normal coordinate
// is negative (meaning that it lies along an edge of the input mesh) and
// returns the halfedge in the input mesh pointing in the same direction
// e.vertex() must live in both meshes
Halfedge
NormalCoordinateIntrinsicTriangulation::getSharedInputEdge(Halfedge he) const {
    verbose_assert(
        intrinsicToInput[he.tailVertex()].type == SurfacePointType::Vertex,
        "I can only identify edges which come out of shared vertices");

    int iE = normalCoordinates.roundabouts[he];
    while (iE < 0) iE += normalCoordinates.roundaboutDegrees[he.vertex()];

    Vertex inputVertex = intrinsicToInput[he.tailVertex()].vertex;
    Halfedge inputHe   = inputVertex.halfedge();

    // step counterclockise iE times
    for (int i = 0; i < iE; i++) inputHe = inputHe.next().next().twin();

    return inputHe;
}

Vertex NormalCoordinateIntrinsicTriangulation::insertVertex(SurfacePoint pt,
                                                            bool verbose) {
    Vertex newVertex;
    switch (pt.type) {
    case SurfacePointType::Vertex:
        newVertex = pt.vertex;
        break;
    case SurfacePointType::Edge:
        edgeSplits++;
        newVertex = splitEdge(pt.edge, pt.tEdge, verbose);
        break;
    case SurfacePointType::Face:
        faceSplits++;
        newVertex = splitFace(pt.face, pt.faceCoords, verbose);
        break;
    }

    return newVertex;
}

Vertex NormalCoordinateIntrinsicTriangulation::splitFace(Face f, Vector3 bary,
                                                         bool verbose) {
    std::clock_t tStart = std::clock();
    // auto data =
    //     normalCoordinates.computeVertexInsertionDataGeodesic(*this, f, bary);
    // Vertex newVertex =
    //     intrinsicMesh->insertVertex(f); // TODO: use mutation Manager
    // cout << "Inserted vertex in face " << f << " at position " << bary <<
    // endl; normalCoordinates.applyVertexInsertionData(newVertex, data); return
    // newVertex;

    std::array<Vector2, 3> vertCoords = vertexCoordinatesInFace(f);
    Vector2 newPCoord = (bary.x * vertCoords[0] + bary.y * vertCoords[1] +
                         bary.z * vertCoords[2]);

    Vector3 fEdgeLengths{
        intrinsicEdgeLengths[f.halfedge().next().edge()],
        intrinsicEdgeLengths[f.halfedge().next().next().edge()],
        intrinsicEdgeLengths[f.halfedge().edge()]};

    std::array<double, 3> newEdgeLengths;
    newEdgeLengths[0] =
        displacementLength(bary - Vector3{1, 0, 0}, fEdgeLengths);
    newEdgeLengths[1] =
        displacementLength(bary - Vector3{0, 1, 0}, fEdgeLengths);
    newEdgeLengths[2] =
        displacementLength(bary - Vector3{0, 0, 1}, fEdgeLengths);
    for (size_t iE = 0; iE < 3; iE++) {
        // newEdgeLengths[iE] = (newPCoord - vertCoords[iE]).norm();
        if (!std::isfinite(newEdgeLengths[iE])) {
            WATCH3(vertCoords, newPCoord, f);
            WATCH(fEdgeLengths);
            throw_verbose_runtime_error("non finite edge length");
        }
        /*
        // Check triangle inequalities
        std::array<double, 3> ol{fEdgeLengths.z, fEdgeLengths.x,
                                 fEdgeLengths.y};
        std::array<double, 3> tl{newEdgeLengths[iE],
                                 newEdgeLengths[(iE + 1) % 3], ol[iE]};
        bool te0 = (tl[0] + tl[1] >= tl[2]);
        bool te1 = (tl[1] + tl[2] >= tl[0]);
        bool te2 = (tl[2] + tl[0] >= tl[1]);

        if (!(te0 && te1 && te2)) {
            std::cerr << "Triangle inequality violation" << std::endl;
            WATCH2(fEdgeLengths, newEdgeLengths);
            WATCH(tl);
            WATCH3(te0, te1, te2);
            return Vertex();
            // throw_verbose_runtime_error(">:(");
        }
        */
    }

    Face insertionFace;
    Vector3 insertionBary;
    std::array<int, 3> counts;

    if (normalCoordinates[f.halfedge().edge()] < 0 &&
        normalCoordinates[f.halfedge().next().edge()] < 0 &&
        normalCoordinates[f.halfedge().next().next().edge()] < 0) {
        // Case 0a: face is empty, all edges shared
        // Note that this means all vertices must be shared
        //

        Halfedge inputHalfedge = identifyInputEdge(f.halfedge());
        insertionFace          = inputHalfedge.face();
        insertionBary          = bary;

        if (inputHalfedge == insertionFace.halfedge()) {
            // good!
        } else if (inputHalfedge == insertionFace.halfedge().next()) {
            insertionBary = rotate(rotate(insertionBary));
        } else if (inputHalfedge == insertionFace.halfedge().next().next()) {
            insertionBary = rotate(insertionBary);
        } else {
            throw_verbose_runtime_error(
                "Face " + std::to_string(insertionFace) + " is not triangular");
        }
        counts = {0, 0, 0};
        easyInsertions1++;
    } else if (normalCoordinates[f.halfedge().edge()] <= 0 &&
               normalCoordinates[f.halfedge().next().edge()] <= 0 &&
               normalCoordinates[f.halfedge().next().next().edge()] <= 0) {
        easyInsertions2++;

        // TODO: use getParentFace
        insertionFace = Face();
        // First, look for a FacePoint - that's the easiest thing to deal with
        for (Vertex v : f.adjacentVertices()) {
            if (intrinsicToInput[v].type == SurfacePointType::Face) {
                if (insertionFace != Face()) {
                    // all faces should agree
                    verbose_assert(insertionFace == intrinsicToInput[v].face,
                                   "if there are no crossings, this face must "
                                   "be contained inside a single input face");
                } else {
                    insertionFace = intrinsicToInput[v].face;
                }
            }
        }

        if (insertionFace == Face()) {
            // If we couldn't find a FacePoint, look for an EdgePoint
            // Note that not all vertices can be shared - if they are, then we
            // can't have any boundary crossings. But this happens sometimes
            // anyway, and is handled next

            auto containsVertex = [](Face f, Vertex v) -> bool {
                for (Vertex vF : f.adjacentVertices()) {
                    if (vF == v) return true;
                }
                return false;
            };

            auto containsEdge = [](Face f, Edge e) -> bool {
                for (Edge eF : f.adjacentEdges()) {
                    if (eF == e) return true;
                }
                return false;
            };

            auto compatible = [&](const SurfacePoint& pt, Face f) -> bool {
                switch (pt.type) {
                case SurfacePointType::Vertex:
                    return containsVertex(f, pt.vertex);
                case SurfacePointType::Edge:
                    return containsEdge(f, pt.edge);
                case SurfacePointType::Face:
                    return pt.face == f;
                }
            };

            for (Vertex v : f.adjacentVertices()) {
                if (intrinsicToInput[v].type == SurfacePointType::Edge) {
                    Edge e  = intrinsicToInput[v].edge;
                    Face f1 = e.halfedge().face();
                    Face f2 = e.halfedge().twin().face();

                    bool f1Okay = e.halfedge().isInterior();
                    bool f2Okay = e.halfedge().twin().isInterior();

                    for (Vertex w : f.adjacentVertices()) {
                        f1Okay = f1Okay && compatible(intrinsicToInput[w], f1);
                        f2Okay = f2Okay && compatible(intrinsicToInput[w], f2);
                    }

                    if ((f1 != f2) && (f1Okay && f2Okay)) {
                        std::cerr
                            << "splitFace err: There are two options for an "
                               "input "
                               "face to insert "
                               "in, and they both look fine. What should I "
                               "do?. (In the meantime, I'm picking option 1)"
                            << myendl;
                    }

                    if (f1Okay) {
                        insertionFace = f1;
                    } else if (f2Okay) {
                        insertionFace = f2;
                    } else {
                        cout << " ==== Failed to find valid face to insert"
                             << myendl;
                        cout << "f1: " << f1
                             << " (isBoundaryLoop: " << f1.isBoundaryLoop()
                             << ")" << myendl;
                        cout << "f2: " << f2
                             << " (isBoundaryLoop: " << f2.isBoundaryLoop()
                             << ")" << myendl;
                        cout << "Underlying edge: " << e
                             << " (isBoundary: " << e.isBoundary() << ")"
                             << myendl;
                        WATCH(e.halfedge().isInterior());
                        WATCH(e.halfedge().twin().isInterior());

                        throw_verbose_runtime_error(
                            "Could not find a valid face in inputMesh to "
                            "insert the new vertex into");
                    }
                }
            }
        }

        if (insertionFace == Face()) {
            // In really degenerate situations, we can have triangles between
            // three
            bool allSharedVertices = true;
            for (Vertex v : f.adjacentVertices()) {
                allSharedVertices =
                    allSharedVertices &&
                    intrinsicToInput[v].type == SurfacePointType::Vertex;
            }
            if (allSharedVertices) {
                // Hope that we find a shared edge. If not, something terrible
                // has happened
                cout << intrinsicToInput[intrinsicMesh->vertex(2336)] << myendl;
                for (Halfedge he : f.adjacentHalfedges()) {
                    WATCH3(he.tailVertex(), normalCoordinates[he.edge()],
                           he.tipVertex());
                    cout << myendl;
                }

                for (Halfedge he : f.adjacentHalfedges()) {
                    if (normalCoordinates[he.edge()] < 0) {
                        Halfedge inputHalfedge = identifyInputEdge(he);
                        insertionFace          = inputHalfedge.face();
                        insertionBary          = bary;

                        if (inputHalfedge == insertionFace.halfedge()) {
                            // good!
                        } else if (inputHalfedge ==
                                   insertionFace.halfedge().next()) {
                            insertionBary = rotate(rotate(insertionBary));
                        } else if (inputHalfedge ==
                                   insertionFace.halfedge().next().next()) {
                            insertionBary = rotate(insertionBary);
                        } else {
                            throw_verbose_runtime_error(
                                "Face " + std::to_string(insertionFace) +
                                " is not triangular");
                        }
                        counts = {0, 0, 0};
                        break;
                    }
                }
            }
        }

        verbose_assert(insertionFace != Face(), "failed to find FacePoint.");

        std::array<Vector3, 3> vertexBary;
        size_t iV = 0;
        for (Vertex v : f.adjacentVertices()) {
            vertexBary[iV] =
                intrinsicToInput[v].inFace(insertionFace).faceCoords;
            iV++;
        }

        insertionBary = bary.x * vertexBary[0] + bary.y * vertexBary[1] +
                        bary.z * vertexBary[2];

        if (verbose) {
            WATCH(vertexBary);
            WATCH(insertionBary);
            for (Vertex v : f.adjacentVertices()) {
                cout << intrinsicToInput[v] << myendl;
                cout << "\t" << intrinsicToInput[v].inFace(insertionFace)
                     << myendl;
            }
        }
        counts = {0, 0, 0};

    } else {
        // Populate the crossing locations for the edges of the triangle
        size_t iHe = 0;
        std::array<std::vector<Curve>, 3> curves;
        std::array<std::vector<std::vector<std::pair<SurfacePoint, double>>>, 3>
            geodesics;
        std::array<std::vector<double>, 3> transverseCrossingTimes;
        std::array<std::vector<double>, 3> boundaryCrossings;
        std::array<std::vector<size_t>, 3> crossingID;
        for (Halfedge he : f.adjacentHalfedges()) {
            for (int ind = 0; ind < normalCoordinates[he.edge()]; ind++) {

                // Get the topological crossings for the curve
                Curve crossings;
                int centerCrossInd;
                std::tie(crossings, centerCrossInd) =
                    normalCoordinates.topologicalTraceBidirectional(he, ind);

                curves[iHe].push_back(crossings);
                geodesics[iHe].push_back(generateFullSingleGeodesicGeometry(
                    *intrinsicMesh, *this, crossings));
                std::vector<std::pair<SurfacePoint, double>>& geodesic =
                    geodesics[iHe].back();
                SurfacePoint& thisCross = geodesic[centerCrossInd + 1].first;
                double tCross           = (he.edge().halfedge() == he)
                                    ? thisCross.tEdge
                                    : 1 - thisCross.tEdge;
                boundaryCrossings[iHe].push_back(tCross);
                transverseCrossingTimes[iHe].push_back(
                    geodesic[centerCrossInd + 1].second);
                crossingID[iHe].push_back(centerCrossInd + 1);
            }
            iHe++;
        }

        hardInsertions++;
        tracingTime += (std::clock() - tStart) / (double)CLOCKS_PER_SEC;
        tStart = std::clock();

        // TODO apply some sanity policies, like that the crossings should be
        // correctly ordered

        counts = computeVertexInsertionCrossingCounts(bary, boundaryCrossings,
                                                      verbose);

        // Find a halfedge bounding the region that our inserted vertex will
        // live inserted in
        std::array<Corner, 3> faceCorners{f.halfedge().corner(),
                                          f.halfedge().next().corner(),
                                          f.halfedge().next().next().corner()};

        auto heIndex = [](Halfedge he) -> int {
            if (he == he.face().halfedge()) {
                return 0;
            } else if (he == he.face().halfedge().next()) {
                return 1;
            } else {
                return 2;
            }
        };

        auto next   = [](int i) { return (i + 1) % 3; };
        auto heBary = [&](Halfedge he, double t) -> Vector3 {
            int i = heIndex(he);
            int j = next(i);

            Vector3 bary = Vector3::zero();
            bary[i]      = (1 - t);
            bary[j]      = (t);

            return bary;
        };

        // Compute the intrinsic and input coordinates of a crossing
        // along a boundary halfedge of intrinsic face f.
        // If pos == 0, returns the intrinsic and input coordinates for the
        // vertex
        auto computeIntrinsicAndInputPoints =
            [&](Halfedge he, int pos) -> std::pair<Vector3, SurfacePoint> {
            if (pos == 0) {
                // source vertex
                return {heBary(he, 0), intrinsicToInput[he.vertex()]};
            } else if (pos == positivePart(normalCoordinates[he.edge()]) + 1) {
                // target vertex
                return {heBary(he, 1), intrinsicToInput[he.next().vertex()]};
            } else {
                // intermediate crossing
                double tInput = transverseCrossingTimes[heIndex(he)][pos - 1];
                Halfedge heInput =
                    identifyInputEdge(curves[heIndex(he)][pos - 1]);
                double tIntrinsic = boundaryCrossings[heIndex(he)][pos - 1];
                return {heBary(he, tIntrinsic), SurfacePoint(heInput, tInput)};
            }
        };

        // Identify 2 points in the input face and their corresponding positions
        // in the intrinsic face
        // TODO: pick 2 points to be spread out. We get numerical errors if the
        // chosen points are too close together

        // Identify as many points as possible on the input face and their
        // corresponding positions on the intrinsic face. In theory 2 suffice,
        // but using more improves stability

        Face inputFace;
        std::vector<std::pair<Vector3, SurfacePoint>> intrinsicInputPairs;
        // Vector3 firstInputPoint, secondInputPoint;
        // Vector3 firstIntrinsicPoint, secondIntrinsicPoint;
        bool foundRegion = false;

        // Used to print error messages later
        size_t myCase = 0;

        // Case 1: in a corner, or just past a corner
        for (size_t iC = 0; iC < 3 && !foundRegion; iC++) {
            int cornerCoord = static_cast<int>(
                normalCoordinates.strictCornerCoord(faceCorners[iC]));

            if (cornerCoord < 1 || counts[iC] > cornerCoord) continue;

            if (verbose) cout << "Case I" << myendl;
            myCase = 1;

            bool justPast =
                (counts[iC] ==
                 static_cast<int>(
                     normalCoordinates.strictCornerCoord(faceCorners[iC])));

            // new vertex is contained in corner iC
            // Bounding halfedge is the counts[iC]'th crossing along
            // corner.halfedge()

            Halfedge firstHedge  = faceCorners[iC].halfedge();
            Halfedge secondHedge = firstHedge.next().next();

            // Grab "outside" halfedge
            int firstInputPointIndex = counts[iC];
            if (justPast) firstInputPointIndex -= 1;
            int secondInputPointIndex = normalCoordinates[secondHedge.edge()] -
                                        (firstInputPointIndex + 1);

            auto& inputHedgeCurve =
                curves[heIndex(firstHedge)][firstInputPointIndex];

            // Identify halfedge in input mesh
            // TODO: repeated effort in computeIntrinsicAndInputPoints
            // Have to take twin due to weird orientation conventions
            Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve).twin();

            if (justPast) { // flip orientation
                inputHalfedge = inputHalfedge.twin();
            }

            inputFace = inputHalfedge.face();

            intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                firstHedge, firstInputPointIndex + 1));
            intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                secondHedge, secondInputPointIndex + 1));

            if (justPast) {
                // If we're just past the end of a corner cell, add in the next
                // crossing along firstHedge and the previous crossing along
                // secondHedge
                // TODO: special case for central polygon?
                intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                    firstHedge, firstInputPointIndex + 2));
                intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                    secondHedge, secondInputPointIndex));
            } else {
                // If we're in a corner cell, add in the previous crossing along
                // firstHedge
                intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                    firstHedge, firstInputPointIndex));

                // If the previous crossing along firstHedge wasn't a vertex,
                // we're in a quad and can also add the next crossing along
                // secondHedge
                if (firstInputPointIndex > 0) {
                    intrinsicInputPairs.push_back(
                        computeIntrinsicAndInputPoints(
                            secondHedge, secondInputPointIndex + 2));
                }
            }

            // for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {
            //     auto p = intrinsicInputPairs[iP];
            //     if (p.second.type == SurfacePointType::Edge &&
            //         p.second.edge.getIndex() == 6295 &&
            //         !checkAdjacent(p.second,
            //                        SurfacePoint(inputFace, Vector3::zero())))
            //                        {
            //         cout << "CASE I" << myendl;
            //         WATCH(iP);
            //         WATCH(intrinsicInputPairs);
            //         WATCH2(firstInputPointIndex, secondInputPointIndex);
            //         WATCH3(normalCoordinates[firstHedge.edge()],
            //                normalCoordinates[firstHedge.next().edge()],
            //                normalCoordinates[secondHedge.edge()]);
            //         break;
            //     }
            // }

            for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {
                if (!checkAdjacent(intrinsicInputPairs[iP].second,
                                   SurfacePoint(inputFace, Vector3::zero()))) {
                    cout << "Case: I, justPast: " << justPast << myendl;
                    cout << "\tfirstInputPointIndex: " << firstInputPointIndex
                         << myendl;
                    cout << "\tsecondInputPointIndex: " << secondInputPointIndex
                         << myendl;
                    cout << "\tpoint: " << intrinsicInputPairs[iP].second
                         << myendl;
                    cout << "\tface: " << inputFace << myendl;
                    cout << "\tiP: " << iP << myendl;
                    for (size_t jP = 0; jP < intrinsicInputPairs.size(); jP++) {
                        cout << "\t\t point " << jP << " : "
                             << intrinsicInputPairs[jP].second << myendl;
                    }
                }
            }

            foundRegion = true;
        }

        // Case 2: in a corner whose corner coordinate is 0
        // Note that we must be in the fan configuration. Otherwise we would
        // also be one-past the final curve crossing some opposite corner
        for (int iC = 0; iC < 3 && !foundRegion; iC++) {
            if (counts[iC] != 0) continue;

            Halfedge firstHedge;
            bool hasFanEdge =
                normalCoordinates.triangleInequalityViolation(f, firstHedge);

            if (!hasFanEdge) {
                WATCH(f);
                for (Edge e : f.adjacentEdges()) {
                    cout << normalCoordinates[e] << myendl;
                }
                WATCH(counts);
                verbose_assert(hasFanEdge,
                               "Triangles in Case 2 must be fan triangles");
            }
            int longHe = heIndex(firstHedge);

            if (iC == next(longHe)) {
                // Covered by case 3
                continue;
            } else if (iC == next(next(longHe))) {
                // next(next(longHe)) is the fan corner, which doesn't work.
                // One of the other corners must be better
                continue;
            }
            if (verbose) cout << "Case II" << myendl;
            myCase = 2;

            // We must have iC == longHe
            // So the new point is in the first section along longHe

            auto& inputHedgeCurve = curves[heIndex(firstHedge)][0];

            // Identify halfedge in input mesh
            Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve).twin();

            inputFace = inputHalfedge.face();

            intrinsicInputPairs.push_back(
                computeIntrinsicAndInputPoints(firstHedge, 1));
            intrinsicInputPairs.push_back(
                computeIntrinsicAndInputPoints(firstHedge, 0));
            intrinsicInputPairs.push_back(
                computeIntrinsicAndInputPoints(firstHedge.next().next(), 0));

            foundRegion = true;
        }

        // Case 3: in the center of the fan region
        if (!foundRegion) {
            if (verbose) cout << "Case III" << myendl;
            myCase = 3;

            Halfedge firstHedge;
            bool hasFanEdge =
                normalCoordinates.triangleInequalityViolation(f, firstHedge);
            verbose_assert(hasFanEdge,
                           "Triangles in Case 3 must be fan triangles");
            int longHe = heIndex(firstHedge);

            // WATCH2(firstHedge, longHe);
            // WATCH(counts);

            // Fan configuration
            int inputPointIndex = counts[longHe] - 1;

            auto& inputHedgeCurve =
                curves[heIndex(firstHedge)][inputPointIndex];

            // Identify halfedge in input mesh
            Halfedge inputHalfedge = identifyInputEdge(inputHedgeCurve);
            inputFace              = inputHalfedge.face();

            intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                firstHedge, inputPointIndex + 1));
            intrinsicInputPairs.push_back(
                computeIntrinsicAndInputPoints(firstHedge.next().next(), 0));

            // Also add in previous crossing
            intrinsicInputPairs.push_back(computeIntrinsicAndInputPoints(
                firstHedge, inputPointIndex + 2));
        }

        // === Compute input position via least squares in barycentric
        // coordinates

        // Use least squares to express intrinsic barycentric coordinate of
        // inserted vertex as a linear combination of barycentric coordinates of
        // intrinsic points on boundary

        // Let P be the matrix if barycentric coordinates for the boundary
        // points
        // We want to find a vector of coefficients b such that P * b = bary,
        // the barycentric coordinate of the input point
        // This is underdetermined, so we find the minimal-norm solution, i.e.
        // min |b|^2 such that Pb = bary

        // Intrinsic barycentric coordinate matrix
        Eigen::MatrixXd P(3, intrinsicInputPairs.size());
        for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {
            Vector3 intBary = intrinsicInputPairs[iP].first;
            P(0, iP)        = intBary.x;
            P(1, iP)        = intBary.y;
            P(2, iP)        = intBary.z;
        }
        Eigen::VectorXd rhs(3);
        rhs << bary.x, bary.y, bary.z;

        Eigen::MatrixXd PPT = P * P.transpose();

        Eigen::VectorXd lambda = PPT.colPivHouseholderQr().solve(rhs);
        Eigen::VectorXd b      = P.transpose() * lambda;

        Eigen::VectorXd err = P * b - rhs;

        insertionBary = Vector3::zero();

        for (size_t iP = 0; iP < intrinsicInputPairs.size(); iP++) {

            // if (!checkAdjacent(intrinsicInputPairs[iP].second,
            //                    SurfacePoint(inputFace, Vector3::zero()))) {
            //     cout << "Case: " << myCase << myendl;
            //     if (myCase == 1) {
            //         cout << "\tjustPast: " << isJustPast << myendl;
            //     }
            // }

            Vector3 inputCrossingBary =
                intrinsicInputPairs[iP].second.inFace(inputFace).faceCoords;
            insertionBary += b(iP) * inputCrossingBary;
        }


        // Very rarely, due to floating point problems, we get barycentric
        // coordinates with one big value and one small value. In this case,
        // we just set them each to be (1-last value)/2
        double roundingTol = 0.5;
        if (insertionBary.x < -roundingTol || insertionBary.y < -roundingTol ||
            insertionBary.z < -roundingTol ||
            insertionBary.x > 1 + roundingTol ||
            insertionBary.y > 1 + roundingTol ||
            insertionBary.z > 1 + roundingTol) {
            WATCH(insertionBary);
            if (insertionBary.x > 1 + roundingTol) {
                if (insertionBary.y < 0) {
                    double s        = (1. - insertionBary.z) / 2.;
                    insertionBary.x = s;
                    insertionBary.y = s;
                } else if (insertionBary.z < 0) {
                    double s        = (1. - insertionBary.y) / 2.;
                    insertionBary.x = s;
                    insertionBary.z = s;
                }
            } else if (insertionBary.y > 1 + roundingTol) {
                if (insertionBary.x < 0) {
                    double s        = (1. - insertionBary.z) / 2.;
                    insertionBary.x = s;
                    insertionBary.y = s;
                } else if (insertionBary.z < 0) {
                    double s        = (1. - insertionBary.x) / 2.;
                    insertionBary.y = s;
                    insertionBary.z = s;
                }
            } else if (insertionBary.z > 1 + roundingTol) {
                if (insertionBary.x < 0) {
                    double s        = (1. - insertionBary.y) / 2.;
                    insertionBary.x = s;
                    insertionBary.z = s;
                } else if (insertionBary.y < 0) {
                    double s        = (1. - insertionBary.x) / 2.;
                    insertionBary.y = s;
                    insertionBary.z = s;
                }
            }
            std::cerr << "insertionBary rounded to " << insertionBary << myendl;
        }

        // HACK: clamp very small values
        // TODO: edge split?
        insertionBary.x = clamp(insertionBary.x, 0., 1.);
        insertionBary.y = clamp(insertionBary.y, 0., 1.);
        insertionBary.z = clamp(insertionBary.z, 0., 1.);
        insertionBary /= (insertionBary.x + insertionBary.y + insertionBary.z);


        insertionFace = inputFace;
    }

    // reorder to fit order of edges incident on newVertex
    std::swap(newEdgeLengths[1], newEdgeLengths[2]);


    auto data = normalCoordinates.computeVertexInsertionData(f, counts);
    // cout << "new2 data: " << data << myendl;

    Vertex newVertex =
        intrinsicMesh->insertVertex(f); // TODO: use mutation Manager
    // cout << "Inserted vertex in face " << f << " at position " << bary
    //      << myendl;
    normalCoordinates.applyVertexInsertionData(newVertex, data);

    SurfacePoint inputPosition(insertionFace, insertionBary);
    intrinsicToInput[newVertex] = inputPosition;

    size_t iE = 0;
    for (Halfedge he : newVertex.outgoingHalfedges()) {
        Edge e                  = he.edge();
        edgeLengths[e]          = newEdgeLengths[iE];
        intrinsicEdgeLengths[e] = newEdgeLengths[iE];

        normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
        normalCoordinates.roundabouts[he] = 0;

        iE++;
    }
    normalCoordinates.roundaboutDegrees[newVertex] = 0;

    // Update quantities
    // refreshQuantities();

    for (Face f : newVertex.adjacentFaces()) {
        // depends on intrinsicEdgeLengths
        updateFaceArea(f);

        // depends on intrinsicEdgeLengths, faceAreas
        updateHalfedgeVectorsInFace(f);

        for (Corner c : f.adjacentCorners()) {
            // depends on intrinsicEdgeLengths
            updateCornerAngle(c);
        }
    }

    for (Vertex v : newVertex.adjacentVertices()) {
        // depends on cornerAngles, vertexAngleSums
        updateHalfedgeVectorsInVertex(v);
    }

    updateVertexAngleSum(newVertex);
    updateHalfedgeVectorsInVertex(newVertex);

    invokeFaceSplitCallbacks(f, newVertex);

    insertionTime += (std::clock() - tStart) / (double)CLOCKS_PER_SEC;
    return newVertex;
}

Vertex NormalCoordinateIntrinsicTriangulation::splitEdge(Edge e, double bary,
                                                         bool verbose) {
    return (e.isBoundary()) ? splitBoundaryEdge(e, bary, verbose)
                            : splitInteriorEdge(e, bary, verbose);
}

Vertex NormalCoordinateIntrinsicTriangulation::splitBoundaryEdge(Edge e,
                                                                 double bary,
                                                                 bool verbose) {
    auto inEdge = [](Edge e, SurfacePoint p) -> SurfacePoint {
        switch (p.type) {
        case SurfacePointType::Vertex:
            if (p.vertex == e.halfedge().tailVertex()) {
                return SurfacePoint(e, 0);
            } else if (p.vertex == e.halfedge().tipVertex()) {
                return SurfacePoint(e, 1);
            }
            break;
        case SurfacePointType::Edge:
            if (p.edge == e) {
                return p;
            }
            break;
        default:
            break;
        }
        throw_verbose_runtime_error("SurfacePoint not in edge");
    };

    auto heBary = [&](Halfedge he, double t) -> Vector3 {
        int i = halfedgeIndexInTriangle(he);
        int j = (i + 1) % 3;

        Vector3 bary = Vector3::zero();
        bary[i]      = (1 - t);
        bary[j]      = (t);

        return bary;
    };
    auto faceEdgeLengths = [&](Face f) -> Vector3 {
        // lengths[i] is the length of the edge opposite the i'th vertex
        return Vector3{intrinsicEdgeLengths[f.halfedge().next().edge()],
                       intrinsicEdgeLengths[f.halfedge().next().next().edge()],
                       intrinsicEdgeLengths[f.halfedge().edge()]};
    };

    if (normalCoordinates[e] >= 0) {
        if (verbose) cout << "Easy Edge Split" << myendl;

        throw_verbose_runtime_error(
            "NSHARP: I think we can get rid of this whole case, since we will "
            "never have normalCoordinates[e] > 0 at boundary. The one case "
            "below for the == 0 case should be sufficient. But it's untested "
            "for == 0, so uncomment this and try it.");

        return Vertex();

        /*
        // Easy case - edge not shared
        // TODO: use normal coordinate edge split code explicitly - it
        // needs fewer geodesic crossing points

        // throw_verbose_runtime_error("Not fully implemented yet - need to mark
        // split segments as fixed if e is");

        Vertex vTipBefore  = e.halfedge().tipVertex();
        Vertex vTailBefore = e.halfedge().tailVertex();
        bool fixedBefore   = isFixed[e];
        isFixed[e]         = false;

        Vertex newVertex =
            splitFace(e.halfedge().face(), heBary(e.halfedge(), bary));
        flipEdgeIfPossible(e);


        // Mark new edges as fixed
        // TODO FIXME this search could potentially fail on a gnarly
        // Delta-complex if the input was a self edge? Fix by passing back the
        // appropriate halfedge from the function below
        Halfedge newHalfedge;
        for (Halfedge he : newVertex.outgoingHalfedges()) {
            if (he.tipVertex() == vTipBefore) {
                newHalfedge = he;
            }
            if (he.tipVertex() == vTipBefore || he.tipVertex() == vTailBefore) {
                if (fixedBefore) {
                    isFixed[he.edge()] = true;
                    // fixedEdges.push_back(e);
                }
            }
        }


        // TODO: invokeEdgeSplitCallbacks here
        // throw_verbose_runtime_error("Didn't call edge split callbacks");

        invokeEdgeSplitCallbacks(
            e, newHalfedge,
            newHalfedge.next().next().twin().next().next().twin());

        return newVertex;
        */

    } else {
        if (verbose) cout << "Shared Edge Split" << myendl;
        // Hard case - edge also exists in input mesh
        Edge inputEdge;
        SurfacePoint src = intrinsicToInput[e.halfedge().tailVertex()];
        SurfacePoint dst = intrinsicToInput[e.halfedge().tipVertex()];

        if (src.type == SurfacePointType::Edge) {
            inputEdge = src.edge;
        } else if (dst.type == SurfacePointType::Edge) {
            inputEdge = dst.edge;
        } else {
            inputEdge = getSharedInputEdge(e.halfedge()).edge();
        }

        double tSrc = inEdge(inputEdge, src).tEdge;
        double tDst = inEdge(inputEdge, dst).tEdge;

        double tInsertion = (1 - bary) * tSrc + bary * tDst;

        SurfacePoint inputPoint(inputEdge, tInsertion);

        std::array<int, 3> newNormalCoordinates =
            normalCoordinates.computeBoundaryEdgeSplitDataGeodesic(*this, e,
                                                                   bary);

        // Compute new edge lengths
        double oldLen = intrinsicEdgeLengths[e];
        std::array<double, 3> newEdgeLengths{(1 - bary) * oldLen, 0,
                                             bary * oldLen};
        newEdgeLengths[1] = displacementLength(
            heBary(e.halfedge(), bary) - heBary(e.halfedge().next(), 1),
            faceEdgeLengths(e.halfedge().face()));

        std::array<bool, 3> newEdgeFixed{isFixed[e], false, isFixed[e]};

        if (verbose) {
            WATCH2(newNormalCoordinates, newEdgeLengths);
        }

        // "edge 2"
        Halfedge newHalfedge =
            intrinsicMesh->splitEdgeTriangular(e); // TODO: use mutation Manager
        Vertex newVertex            = newHalfedge.vertex();
        intrinsicToInput[newVertex] = inputPoint;

        verbose_assert(newHalfedge.isInterior() &&
                           !newHalfedge.twin().isInterior(),
                       "I'm wrong about orientation conventions");

        if (verbose) {
            WATCH(inputPoint);
        }

        // TODO: write applyVertexInsertionData in NormalCoordinates
        // class

        size_t iE = 0; // TODO: check indexing convention
        // Explicit loop to go counterclockwise
        Halfedge he = newHalfedge;
        while (true) {
            Edge e                          = he.edge();
            edgeLengths[e]                  = newEdgeLengths[iE];
            intrinsicEdgeLengths[e]         = newEdgeLengths[iE];
            normalCoordinates.edgeCoords[e] = newNormalCoordinates[iE];

            if (!isFixed[e] && newEdgeFixed[iE]) {
                isFixed[e] = true;
                // fixedEdges.push_back(e);
            } else if (isFixed[e] && !newEdgeFixed[iE]) {
                throw_verbose_runtime_error(
                    "Need to remove e from the fixed list");
            }

            if (!he.isInterior()) break;
            he = he.next().next().twin();
            iE++;
        }

        for (Halfedge he : newVertex.outgoingHalfedges()) {
            // depends on normalCoordinates.edgeCoords
            normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
            normalCoordinates.roundabouts[he] = 0;
        }
        normalCoordinates.roundaboutDegrees[newVertex] = 0;

        // Update quantities

        for (Face f : newVertex.adjacentFaces()) {
            // depends on intrinsicEdgeLengths
            updateFaceArea(f);

            // depends on intrinsicEdgeLengths, faceAreas
            updateHalfedgeVectorsInFace(f);

            for (Corner c : f.adjacentCorners()) {
                // depends on intrinsicEdgeLengths
                updateCornerAngle(c);
            }
        }

        for (Vertex v : newVertex.adjacentVertices()) {
            // depends on cornerAngles, vertexAngleSums
            updateHalfedgeVectorsInVertex(v);
        }

        updateVertexAngleSum(newVertex);
        updateHalfedgeVectorsInVertex(newVertex);

        invokeEdgeSplitCallbacks(
            e, newHalfedge,
            newHalfedge.next().next().twin().next().next().twin());

        return newVertex;
    }
}

Vertex NormalCoordinateIntrinsicTriangulation::splitInteriorEdge(Edge e,
                                                                 double bary,
                                                                 bool verbose) {
    auto inEdge = [](Edge e, SurfacePoint p) -> SurfacePoint {
        switch (p.type) {
        case SurfacePointType::Vertex:
            if (p.vertex == e.halfedge().tailVertex()) {
                return SurfacePoint(e, 0);
            } else if (p.vertex == e.halfedge().tipVertex()) {
                return SurfacePoint(e, 1);
            }
            break;
        case SurfacePointType::Edge:
            if (p.edge == e) {
                return p;
            }
            break;
        default:
            break;
        }
        throw_verbose_runtime_error("SurfacePoint not in edge");
    };

    auto heBary = [&](Halfedge he, double t) -> Vector3 {
        int i = halfedgeIndexInTriangle(he);
        int j = (i + 1) % 3;

        Vector3 bary = Vector3::zero();
        bary[i]      = (1 - t);
        bary[j]      = (t);

        return bary;
    };
    auto faceEdgeLengths = [&](Face f) -> Vector3 {
        // lengths[i] is the length of the edge opposite the i'th vertex
        return Vector3{intrinsicEdgeLengths[f.halfedge().next().edge()],
                       intrinsicEdgeLengths[f.halfedge().next().next().edge()],
                       intrinsicEdgeLengths[f.halfedge().edge()]};
    };

    if (normalCoordinates[e] >= 0) {
        if (verbose) cout << "Easy Edge Split" << myendl;
        // Easy case - edge not shared
        // TODO: use normal coordinate edge split code explicitly - it
        // needs fewer geodesic crossing points

        Vertex vTipBefore  = e.halfedge().tipVertex();
        Vertex vTailBefore = e.halfedge().tailVertex();
        bool fixedBefore   = isFixed[e];

        Vertex newVertex =
            splitFace(e.halfedge().face(), heBary(e.halfedge(), bary));
        bool flipHappened = flipEdgeIfPossible(e);

        if (!flipHappened)
            throw_verbose_runtime_error("pos-split flip failed!");

        // Find two of the halfedges that together make up the old edge
        Halfedge newHalfedge = e.halfedge().prevOrbitFace().twin();
        Halfedge otherNewHalfedge =
            newHalfedge.next().next().twin().next().next().twin();


        // Update is-fixed array
        if (fixedBefore) {
            isFixed[newHalfedge.edge()]      = true;
            isFixed[otherNewHalfedge.edge()] = true;
        }

        invokeEdgeSplitCallbacks(e, newHalfedge, otherNewHalfedge);

        return newVertex;

    } else {
        if (verbose) cout << "Shared Edge Split" << myendl;
        // Hard case - edge also exists in input mesh
        Edge inputEdge;
        SurfacePoint src = intrinsicToInput[e.halfedge().tailVertex()];
        SurfacePoint dst = intrinsicToInput[e.halfedge().tipVertex()];

        if (src.type == SurfacePointType::Edge) {
            inputEdge = src.edge;
        } else if (dst.type == SurfacePointType::Edge) {
            inputEdge = dst.edge;
        } else {
            inputEdge = getSharedInputEdge(e.halfedge()).edge();
        }

        double tSrc = inEdge(inputEdge, src).tEdge;
        double tDst = inEdge(inputEdge, dst).tEdge;

        double tInsertion = (1 - bary) * tSrc + bary * tDst;

        SurfacePoint inputPoint(inputEdge, tInsertion);

        std::array<int, 4> newNormalCoordinates =
            normalCoordinates.computeInteriorEdgeSplitDataGeodesic(*this, e,
                                                                   bary);

        // Compute new edge lengths
        double oldLen = intrinsicEdgeLengths[e];
        std::array<double, 4> newEdgeLengths{0, (1 - bary) * oldLen, 0,
                                             bary * oldLen};
        newEdgeLengths[0] =
            displacementLength(heBary(e.halfedge().twin(), 1 - bary) -
                                   heBary(e.halfedge().twin().next(), 1),
                               faceEdgeLengths(e.halfedge().twin().face()));
        newEdgeLengths[2] = displacementLength(
            heBary(e.halfedge(), bary) - heBary(e.halfedge().next(), 1),
            faceEdgeLengths(e.halfedge().face()));

        std::array<bool, 4> newEdgeFixed{false, isFixed[e], false, isFixed[e]};

        if (verbose) {
            WATCH2(newNormalCoordinates, newEdgeLengths);
        }

        // "edge 2"
        Halfedge newHalfedge =
            intrinsicMesh->splitEdgeTriangular(e); // TODO: use mutation Manager
        Vertex newVertex            = newHalfedge.vertex();
        intrinsicToInput[newVertex] = inputPoint;

        if (verbose) {
            WATCH(inputPoint);
        }

        // TODO: write applyVertexInsertionData in NormalCoordinates
        // class

        size_t iE = 1; // TODO: fix indexing convention
        for (Halfedge he : newVertex.outgoingHalfedges()) {
            Edge e = he.edge();

            edgeLengths[e]                  = newEdgeLengths[iE];
            intrinsicEdgeLengths[e]         = newEdgeLengths[iE];
            normalCoordinates.edgeCoords[e] = newNormalCoordinates[iE];

            if (!isFixed[e] && newEdgeFixed[iE]) {
                isFixed[e] = true;
                // fixedEdges.push_back(e);
            } else if (isFixed[e] && !newEdgeFixed[iE]) {
                throw_verbose_runtime_error(
                    "Need to remove e from the fixed list");
            }

            // indexing goes counterclockwise, but loop goes clockwise
            iE = (iE + 3) % 4;
        }

        for (Halfedge he : newVertex.outgoingHalfedges()) {
            // depends on normalCoordinates.edgeCoords
            normalCoordinates.setRoundaboutFromPrevRoundabout(he.twin());
            normalCoordinates.roundabouts[he] = 0;
        }
        normalCoordinates.roundaboutDegrees[newVertex] = 0;

        // Update quantities

        for (Face f : newVertex.adjacentFaces()) {
            // depends on intrinsicEdgeLengths
            updateFaceArea(f);

            // depends on intrinsicEdgeLengths, faceAreas
            updateHalfedgeVectorsInFace(f);

            for (Corner c : f.adjacentCorners()) {
                // depends on intrinsicEdgeLengths
                updateCornerAngle(c);
            }
        }

        for (Vertex v : newVertex.adjacentVertices()) {
            // depends on cornerAngles, vertexAngleSums
            updateHalfedgeVectorsInVertex(v);
        }

        updateVertexAngleSum(newVertex);
        updateHalfedgeVectorsInVertex(newVertex);

        invokeEdgeSplitCallbacks(
            e, newHalfedge,
            newHalfedge.next().next().twin().next().next().twin());

        return newVertex;
    }
}

Face NormalCoordinateIntrinsicTriangulation::removeInsertedVertex(Vertex v) {
    // Stolen from geometrycentral/signpost_intrinsic_triangulation.cpp
    // Strategy: flip edges until the vertex has degree three, then remove by
    // replacing with a single face
    // TODO needs a proof that this always works... what about self edges, etc?
    // Seems to work well.

    // What about starting with degree < 3? Since this vertex necessarily has
    // angle sum 2PI, this could only happen in the case of degree 2, with
    // exactly degenerate triangles. Since we assume non-degenerate triangles
    // throughout, we'll consider that to not happen.

    if (intrinsicToInput[v].type == SurfacePointType::Vertex)
        return Face(); // can't remove original vertices

    if (isOnFixedEdge(v)) {
        return Face(); // don't try to remove boundary vertices, for now at
                       // least
    }

    if (v.isBoundary()) {
        throw_verbose_runtime_error("boundary vertex removal not implemented");
    }

    // === Flip edges until v has degree 3


    size_t iterCount = 0;
    while (v.degree() != 3 && iterCount < 10 * v.degree()) {

        // Find the highest priority edge to flip
        Edge bestFlipEdge;
        double bestFlipScore = -std::numeric_limits<double>::infinity();
        bool bestFlipIsLoop  = false;
        for (Edge e : v.adjacentEdges()) {

            double flipScore = checkFlip(e);
            bool isLoop      = e.firstVertex() == e.secondVertex();

            // This logic picks the most-preferred edge to flip. The policy is
            // basically "pick the edge whith the highest flipScore", except
            // that we prefer loop edges if there are any.
            if (isLoop) {
                if (bestFlipIsLoop) {
                    // if the one we currently have is a loop, only take this
                    // one if it is better
                    if (flipScore > bestFlipScore) {
                        bestFlipScore = flipScore;
                        bestFlipEdge  = e;
                    }

                } else {
                    // if the one we currently have is not a loop, always take
                    // this one if it is valid
                    if (flipScore > 0.) {
                        bestFlipScore = flipScore;
                        bestFlipEdge  = e;
                    }
                }

                bestFlipIsLoop = true;
            } else {
                if (!bestFlipIsLoop) { // only overwrite if the best is not a
                                       // loop
                    if (flipScore > bestFlipScore) {
                        bestFlipScore  = flipScore;
                        bestFlipEdge   = e;
                        bestFlipIsLoop = false;
                    }
                }
            }
        }

        if (bestFlipEdge == Edge()) {
            return Face();
            // throw_verbose_runtime_error("failed to remove vertex " +
            //                             std::to_string(v) +
            //                             ".  Could not find any edge to
            //                             flip");
        }

        // Passing -inf as the tolerance forces us to always do the flip, since
        // we've already verified it above
        // flipEdgeIfPossible(bestFlipEdge,
        //                    -std::numeric_limits<double>::infinity());
        flipEdgeIfPossible(bestFlipEdge);

        iterCount++;
    }

    // give up if something went wrong (eg. flipped edges)
    if (v.degree() != 3) {
        // throw_verbose_runtime_error(
        //     "failed to remove vertex " + std::to_string(v) +
        //     ".  Somehow vertex degree is not 3. Was it 2 beforehand? (which "
        //     "can only be degenerate)");
        return Face();
    }

    // ==== Remove the vertex
    Face newF = intrinsicMesh->removeVertex(v);

    // ==== Update cached data
    // Edge lengths, normal coordinates, and roundabouts should be okay

    for (Corner c : newF.adjacentCorners()) {
        updateCornerAngle(c);
    }
    updateFaceArea(newF);
    updateHalfedgeVectorsInFace(newF);
    for (Halfedge he : newF.adjacentHalfedges()) {
        updateHalfedgeVectorsInVertex(he.vertex());
    }

    return newF;
}


Vertex NormalCoordinateIntrinsicTriangulation::moveVertex(Vertex v,
                                                          Vector2 vec) {

    // Find the insertion location
    TraceOptions options;
    SurfacePoint startP(v);
    TraceGeodesicResult intrinsicTraceResult =
        traceGeodesic(*this, startP, vec, options);
    SurfacePoint newPositionOnIntrinsic = intrinsicTraceResult.endPoint;
    if (newPositionOnIntrinsic.type == SurfacePointType::Vertex &&
        newPositionOnIntrinsic.vertex == Vertex()) {
        // tracing failed
        return Vertex();
    }

    // Insert a new vertex
    Vertex newV = insertVertex(newPositionOnIntrinsic);

    // Delete the old vertex
    removeInsertedVertex(v);

    return newV;
}

bool NormalCoordinateIntrinsicTriangulation::isDelaunay(Edge e,
                                                        double tol) const {
    double cWeight = getEdgeCotanWeight(e);
    return e.isBoundary() || !(cWeight < -tol);
}

size_t NormalCoordinateIntrinsicTriangulation::nSubdividedVertices() const {
    size_t nNewVertices = 0;
    for (Edge e : intrinsicMesh->edges())
        nNewVertices += positivePart(normalCoordinates[e]);

    return intrinsicMesh->nVertices() + nNewVertices;
}

// HACK: represents arcs parallel to a mesh edge with a single pair {-n,
// he} where n is the number of arcs parallel to he.edge() Trace an edge
// of the input mesh over the intrinsic triangulation
CompoundCurve
NormalCoordinateIntrinsicTriangulation::traceInputEdge(Edge e,
                                                       bool verbose) const {

    // verbose = verbose || e.getIndex() == 18027;

    auto vertexHalfedge = [&](Vertex v, size_t iH) {
        Halfedge he = v.halfedge();

        // Iterate counterclockwise
        for (size_t i = 0; i < iH; ++i) {
            he = he.next().next().twin();
        }
        return he;
    };

    auto traceFollowingComponents =
        [&](const Curve& firstComponent) -> CompoundCurve {
        CompoundCurve cc;
        bool valid          = true;
        Curve nextComponent = firstComponent;
        if (verbose) {
            WATCH(firstComponent);
        }
        while (valid) {
            cc.components.push_back(nextComponent);
            std::tie(valid, nextComponent) =
                traceNextCurve(cc.components.back(), verbose);
            if (verbose) {
                WATCH2(valid, nextComponent);
            }
        }

        return cc;
    };

    auto directPath = [&](Halfedge he) {
        return Curve{{std::make_pair(-1, he)}};
    };

    Vertex vTrace = src(e);

    // TODO: Allow different vertex sets
    Vertex v = intrinsicMesh->vertex(vTrace.getIndex());

    Halfedge he = v.halfedge();


    // Loop over all halfedges of intrinsicMesh coming out of v until we
    // find the one whose corner contains the halfedge e.halfedge() of
    // inputMesh
    // TODO: this could be optimized. Calling vertexHalfedge repeatedly
    // repeats a lot of work
    do {
        size_t iStart = normalCoordinates.roundabouts[he];

        size_t em       = normalCoordinates.strictDegree(he.corner());
        size_t startInd = -negativePart(normalCoordinates[he.edge()]);
        size_t endInd =
            -negativePart(normalCoordinates[he.next().next().edge()]);

        size_t width = em + startInd + endInd;

        // Loop over all halfedges of inputMesh coming out of this
        // corner
        for (size_t iH = 0; iH < width; ++iH) {
            Halfedge heTrace = vertexHalfedge(vTrace, iStart + iH);
            if (heTrace != e.halfedge()) continue;

            if (verbose) {
                cout << "Tracing from vertex " << v << myendl;
                cout << "Found edge " << e << myendl;
                cout << "\t iH = " << iH << " of " << width << " = " << startInd
                     << " + " << em << " + " << endInd << myendl;
                cout << "\t\t index is " << iStart + iH
                     << " from source vertex " << vTrace << myendl;
                cout << "\t\t past halfedge  " << he << " / edge " << he.edge()
                     << myendl;
                cout << "\t ALL HEDGES: " << myendl;
                Halfedge he = v.halfedge();

                // Iterate counterclockwise
                for (size_t i = 0; i < v.degree(); ++i) {
                    cout << "\t\t " << he << "\t (" << he.edge()
                         << ") : " << normalCoordinates.roundabouts[he]
                         << " | edge coord: "
                         << normalCoordinates.edgeCoords[he.edge()]
                         << " | corner deg: "
                         << normalCoordinates.strictDegree(he.corner())
                         << myendl;

                    he = he.next().next().twin();
                }
            }

            if (iH >= width - endInd) {
                Halfedge pathHe = he.next().next().twin();
                return traceFollowingComponents(directPath(pathHe));
            } else if (iH < startInd) {
                return traceFollowingComponents(directPath(he));
            } else {
                size_t idx = iH - startInd;

                return traceFollowingComponents(
                    normalCoordinates.topologicalTrace(
                        he.next(),
                        positivePart(normalCoordinates[he.edge()]) + idx));
            }
        }
        // orbit counterclockwise
        he = he.next().next().twin();
    } while (he != v.halfedge());


    std::cerr << "Something somewhere went horribly wrong" << std::endl;
    he = v.halfedge();
    do {
        size_t iStart = normalCoordinates.roundabouts[he];

        size_t em       = positivePart(normalCoordinates[he.next().edge()] -
                                 normalCoordinates[he.edge()] -
                                 normalCoordinates[he.next().next().edge()]);
        size_t startInd = (normalCoordinates[he.edge()] == 0) ? 1 : 0;
        size_t endInd =
            (normalCoordinates[he.next().next().edge()] == 0) ? 1 : 0;

        size_t width = em + startInd + endInd;
        std::cerr << "iStart: " << iStart << "\tem: " << em
                  << "\tstartInd: " << startInd << "\t endInd: " << endInd
                  << myendl;
        for (size_t iH = 0; iH < width; ++iH) {
            std::cerr << "\ttrying halfedge " << iStart + iH << " of "
                      << vTrace.degree() << " (aka "
                      << normalCoordinates.roundaboutDegrees[v] << ")"
                      << myendl;
        }
        // orbit counterclockwise
        he = he.next().next().twin();
    } while (he != v.halfedge());

    return {};
}

std::pair<bool, Curve>
NormalCoordinateIntrinsicTriangulation::traceNextCurve(const Curve& oldCurve,
                                                       bool verbose) const {
    verbose_assert(!oldCurve.crossings.empty(),
                   "curves must contain some crossings. Even shared edges "
                   "contain a 'parallel' crossing");
    Halfedge finalHe;
    int finalPos;
    std::tie(finalPos, finalHe) = oldCurve.crossings.back();

    if (finalPos < 0) {
        if (verbose) {
            cout << "SHARED EDGE" << myendl;
            WATCH3(finalHe, intrinsicToInput[finalHe.tailVertex()],
                   intrinsicToInput[finalHe.tipVertex()]);
        }

        // Shared edge
        Vertex v              = finalHe.tipVertex();
        SurfacePoint inputPos = intrinsicToInput[v];
        if (inputPos.type == SurfacePointType::Edge) {
            // keep tracing

            // Look for another shared edge
            for (Halfedge he : v.outgoingHalfedges()) {
                if (he != finalHe.twin() && normalCoordinates[he.edge()] < 0) {
                    return {true, Curve{{{-1, he}}}};
                }
            }

            // Look for a non-shared curve
            for (Corner c : v.adjacentCorners()) {
                if (normalCoordinates.strictDegree(c) > 0) {
                    verbose_assert(normalCoordinates.strictDegree(c) == 1,
                                   "there can only be one");
                    // we found a curve. And it's different than
                    // oldCurve since it's not a shared edge
                    return {true, normalCoordinates.topologicalTrace(c, 0)};
                }
            }
        } else {
            return {false, Curve{}};
        }
    } else {
        if (verbose) cout << "Normal Edge" << myendl;
        // Ordinary (transverse) edge
        Corner incomingCorner = finalHe.twin().next().next().corner();
        Vertex v              = incomingCorner.vertex();
        SurfacePoint inputPos = intrinsicToInput[v];
        if (inputPos.type == SurfacePointType::Edge) {
            // keep tracing

            // Look for a shared edge
            for (Halfedge he : v.outgoingHalfedges()) {
                if (normalCoordinates[he.edge()] < 0) {
                    return {true, Curve{{{-1, he}}}};
                }
            }

            // Look for a non-shared curve
            for (Corner c : v.adjacentCorners()) {
                if (normalCoordinates.strictDegree(c) > 0 &&
                    c != incomingCorner) {
                    verbose_assert(normalCoordinates.strictDegree(c) == 1,
                                   "there can only be one");
                    // we found a curve. And it's different than
                    // oldCurve since it's not a shared edge
                    return {true, normalCoordinates.topologicalTrace(c, 0)};
                } else if (normalCoordinates.strictDegree(c) == 2) {
                    throw_verbose_runtime_error(
                        "Your geometry is really really bad");
                    verbose_assert(c == incomingCorner,
                                   "there can only be one");

                    return {true, normalCoordinates.topologicalTrace(
                                      c.halfedge().next(), finalPos)};

                } else {
                    throw_verbose_runtime_error("this shouldn't be possible");
                    return {false, Curve{}};
                }
            }
        } else {
            return {false, Curve{}};
        }
    }
    throw_verbose_runtime_error("This should be unreachable");
}

// Inverse to traceInputEdge
Halfedge
NormalCoordinateIntrinsicTriangulation::identifyInputEdge(const Curve& path,
                                                          bool verbose) const {
    verbose_assert(path.crossings.size() >= 1,
                   "paths of length 0 do not correspond to edges");

    auto vertexHalfedge = [&](Vertex v, size_t iH) {
        Halfedge he = v.halfedge();

        // Iterate counterclockwise
        for (size_t i = 0; i < iH; ++i) {
            he = he.next().next().twin();
        }
        return he;
    };

    if (path.crossings[0].first < 0) {
        // If edge is shared, use roundabouts to identify shared edge

        Halfedge he     = path.crossings[0].second;
        Vertex inputSrc = intrinsicToInput[he.vertex()].vertex;
        return vertexHalfedge(inputSrc, normalCoordinates.roundabouts[he]);
    } else {
        Halfedge he = path.crossings[0].second;
        int p       = path.crossings[0].first -
                normalCoordinates.strictCornerCoord(he.corner());

        if (verbose) {
            std::cerr << "halfedge: " << he
                      << "\t crossing No: " << path.crossings[0].first
                      << "\t cornerCoord: "
                      << normalCoordinates.strictCornerCoord(he.corner())
                      << "\tp: " << p << std::endl;
        }
        verbose_assert(p >= 0, "crossing must have nonnegative index");

        SurfacePoint src = intrinsicToInput[he.next().next().vertex()];
        verbose_assert(src.type == SurfacePointType::Vertex,
                       "edge must start at vertex");
        Vertex inputSrc = src.vertex;

        Halfedge hePrev = he.next().next();
        int rIdx        = p + normalCoordinates.roundabouts[hePrev] -
                   negativePart(normalCoordinates[hePrev.edge()]);

        return vertexHalfedge(inputSrc, rIdx);
    }
}

// Identify shared halfedge, throw exception if halfedge is not shared
// (i.e. edgeCoords[he.edge()] must be negative)
Halfedge
NormalCoordinateIntrinsicTriangulation::identifyInputEdge(Halfedge he) const {
    verbose_assert(normalCoordinates[he.edge()] < 0,
                   "shared edge must have edgeCoord -1");

    auto vertexHalfedge = [&](Vertex v, size_t iH) {
        Halfedge he = v.halfedge();

        // Iterate counterclockwise
        for (size_t i = 0; i < iH; ++i) {
            he = he.next().next().twin();
        }
        return he;
    };

    Vertex src = intrinsicToInput[he.vertex()].vertex;
    return vertexHalfedge(src, normalCoordinates.roundabouts[he]);
}

std::array<Vector2, 3>
NormalCoordinateIntrinsicTriangulation::vertexCoordinatesInFace(
    Face face) const {
    // stolen from
    // gc/intrinsic_geometry_interface.cpp:computeHalfedgeVectorsInFace

    // Gather some values
    Halfedge heAB = face.halfedge();
    Halfedge heBC = heAB.next();
    Halfedge heCA = heBC.next();
    GC_SAFETY_ASSERT(heCA.next() == heAB, "faces must be triangular");

    double lAB = intrinsicEdgeLengths[heAB.edge()];
    double lBC = intrinsicEdgeLengths[heBC.edge()];
    double lCA = intrinsicEdgeLengths[heCA.edge()];

    // Assign positions to all three vertices
    Vector2 pA{0., 0.}; // used implicitly
    Vector2 pB{lAB, 0.};
    // pC is the hard one:

    // Herons formula
    // stolen from
    // gc/intrinsic_geometry_interface.cpp:computeFaceAreas
    double s     = (lAB + lBC + lCA) / 2.0;
    double arg   = s * (s - lAB) * (s - lBC) * (s - lCA);
    arg          = std::fmax(0., arg);
    double tArea = std::sqrt(arg);

    // Compute width and height of right triangle formed via altitude
    // from C
    double h = 2. * tArea / lAB;
    double w = std::sqrt(std::max(0., lCA * lCA - h * h));

    // Take the closer of the positive and negative solutions
    if (lBC * lBC > (lAB * lAB + lCA * lCA)) w *= -1.0;

    // Project some vectors to get the actual position
    Vector2 pC{w, h};

    return {pA, pB, pC};
}

void NormalCoordinateIntrinsicTriangulation::setFixedEdges(
    const EdgeData<bool>& fixedEdges_) {

    isFixed = fixedEdges_;


    // fixedEdges.clear();
    // for (Edge e : intrinsicMesh->edges()) {
    // if (isFixed[e]) fixedEdges.push_back(e);
    //}
}

// If f is entirely contained in some face of the input mesh, return that
// face Otherwise return Face()
Face NormalCoordinateIntrinsicTriangulation::getParentFace(Face f) const {
    auto containsVertex = [](Face f, Vertex v) -> bool {
        for (Vertex vF : f.adjacentVertices()) {
            if (vF == v) return true;
        }
        return false;
    };

    auto containsEdge = [](Face f, Edge e) -> bool {
        for (Edge eF : f.adjacentEdges()) {
            if (eF == e) return true;
        }
        return false;
    };

    auto compatible = [&](const SurfacePoint& pt, Face f) -> bool {
        switch (pt.type) {
        case SurfacePointType::Vertex:
            return containsVertex(f, pt.vertex);
        case SurfacePointType::Edge:
            return containsEdge(f, pt.edge);
        case SurfacePointType::Face:
            return pt.face == f;
        }
    };

    // Look for a FacePoint
    for (Vertex v : f.adjacentVertices()) {
        SurfacePoint vP = intrinsicToInput[v];
        if (vP.type == SurfacePointType::Face) {
            Face parentFace = vP.face;

            // Check if this works for everyone else
            for (Vertex w : f.adjacentVertices()) {
                if (!compatible(intrinsicToInput[w], parentFace)) {
                    return Face();
                }
            }

            return parentFace;
        }
    }

    // Look for an EdgePoint
    for (Vertex v : f.adjacentVertices()) {
        if (intrinsicToInput[v].type == SurfacePointType::Edge) {
            Edge e  = intrinsicToInput[v].edge;
            Face f1 = e.halfedge().face();
            Face f2 = e.halfedge().twin().face();

            bool f1Okay = e.halfedge().isInterior();
            bool f2Okay = e.halfedge().twin().isInterior();

            for (Vertex w : f.adjacentVertices()) {
                f1Okay = f1Okay && compatible(intrinsicToInput[w], f1);
                f2Okay = f2Okay && compatible(intrinsicToInput[w], f2);
            }

            return (f1Okay) ? f1 : (f2Okay) ? f2 : Face();
        }
    }

    // Give up
    return Face();
}

// ======================================================
//          Geometry and Helpers
// ======================================================

void NormalCoordinateIntrinsicTriangulation::computeEdgeLengths() {
    edgeLengths = intrinsicEdgeLengths;
}

void NormalCoordinateIntrinsicTriangulation::invokeEdgeFlipCallbacks(Edge e) {
    for (auto& fn : edgeFlipCallbackList) {
        fn(e);
    }
}
void NormalCoordinateIntrinsicTriangulation::invokeFaceSplitCallbacks(
    Face f, Vertex v) {
    for (auto& fn : faceSplitCallbackList) {
        fn(f, v);
    }
}
void NormalCoordinateIntrinsicTriangulation::invokeEdgeSplitCallbacks(
    Edge e, Halfedge he1, Halfedge he2) {
    for (auto& fn : edgeSplitCallbackList) {
        fn(e, he1, he2);
    }
}

// Compute the cotan weight of halfedge ij in terms of the lengths of
// its neighbors
double halfedgeCotanWeight(double l_ij, double l_jk, double l_ki) {
    double areaV    = triangleArea(l_ij, l_jk, l_ki);
    double cotValue = (-l_ij * l_ij + l_jk * l_jk + l_ki * l_ki) / (4. * areaV);
    return cotValue / 2;
}

FaceData<Vector2>
interpolateTangentVectorsB(const NormalCoordinateIntrinsicTriangulation& tri,
                           const CommonSubdivision& cs,
                           const FaceData<Vector2>& dataB) {

    FaceData<Vector2> interp(*cs.mesh);

    for (Face f : cs.mesh->faces()) {
        Face fB = cs.sourceFaceB[f];
        if (fB == Face()) {
            std::cerr << "Encountered Face() as parent?" << myendl;
            continue;
        }

        // Find the position's of fB's vertices in its tangent space
        // Vector2 p0{0, 0};
        Vector2 p1 = tri.halfedgeVectorsInFace[fB.halfedge()];
        Vector2 p2 = -tri.halfedgeVectorsInFace[fB.halfedge().next().next()];

        // Find the position's of f's vertices in barycentric coords on fB
        Vertex v0      = f.halfedge().tailVertex();
        Vector3 v0Bary = cs.sourcePoints[v0]->posB.inFace(fB).faceCoords;
        Vertex v1      = f.halfedge().tipVertex();
        Vector3 v1Bary = cs.sourcePoints[v1]->posB.inFace(fB).faceCoords;

        // Find the position's of f's vertices in fB's tangent space
        Vector2 q0 = v0Bary.y * p1 + v0Bary.z * p2;
        Vector2 q1 = v1Bary.y * p1 + v1Bary.z * p2;

        //
        Vector2 fBasisInFB = (q1 - q0).normalize();

        interp[f] = dataB[fB] / fBasisInFB;
    }

    return interp;
}
