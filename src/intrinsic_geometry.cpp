#include "geometrycentral/intrinsic_geometry.h"

#include "geometrycentral/discrete_operators.h"

#include <fstream>
#include <limits>

namespace geometrycentral {

 IntrinsicGeometry::IntrinsicGeometry(HalfedgeMesh* mesh_) : mesh(mesh_) {}

 void IntrinsicGeometry::buildDependencies() {

  // Helper to add a quantity, binding this instance to its compute function and adding it to the list of all quantities
  auto addQuantity = [&](DependentQuantity& q, std::vector<DependentQuantity*> deps,
                         std::function< void(IntrinsicGeometry*)> func) {
    q = DependentQuantity(deps, std::bind(func, this));
    allQuantities.push_back(&q);
  };


  // === ALL the quantities
  // clang-format off
  //          quantity manager            dependencies                                  compute function
  //
  // == Basic geometric quantities
  addQuantity(faceAreasQ,                 {&edgeLengthsQ},                              &IntrinsicGeometry::computeFaceAreas);
  addQuantity(vertexDualAreasQ,           {&faceAreasQ},                                &IntrinsicGeometry::computeVertexDualAreas);
  addQuantity(edgeLengthsQ,               {},                                           &IntrinsicGeometry::computeEdgeLengths);
  addQuantity(halfedgeCotanWeightsQ,      {&halfedgeOppositeAnglesQ},                   &IntrinsicGeometry::computeHalfedgeCotanWeights);
  addQuantity(edgeCotanWeightsQ,          {&halfedgeCotanWeightsQ},                     &IntrinsicGeometry::computeEdgeCotanWeights);
  addQuantity(vertexAngleDefectsQ,        {&halfedgeOppositeAnglesQ},                   &IntrinsicGeometry::computeVertexAngleDefects);
  
  // == Vector fields, angles, and transport
  addQuantity(halfedgeFaceCoordsQ,              {&halfedgeOppositeAnglesQ, &edgeLengthsQ},          &IntrinsicGeometry::computeHalfedgeFaceCoords);
  addQuantity(faceTransportCoefsQ,              {&halfedgeFaceCoordsQ},                             &IntrinsicGeometry::computeFaceTransportCoefs);
  addQuantity(halfedgeOppositeAnglesQ,          {&edgeLengthsQ},                                    &IntrinsicGeometry::computeHalfedgeOppositeAngles);
  addQuantity(halfedgeRescaledOppositeAnglesQ,  {&vertexAngleDefectsQ, &halfedgeOppositeAnglesQ},   &IntrinsicGeometry::computeHalfedgeRescaledOppositeAngles);
  addQuantity(halfedgeVertexCoordsQ,            {&halfedgeRescaledOppositeAnglesQ},                 &IntrinsicGeometry::computeHalfedgeVertexCoords);
  addQuantity(vertexTransportCoefsQ,            {&halfedgeVertexCoordsQ},                           &IntrinsicGeometry::computeVertexTransportCoefs);
  
  // == Indices
  addQuantity(vertexIndicesQ,             {},                                           &IntrinsicGeometry::computeVertexIndices);
  addQuantity(interiorVertexIndicesQ,     {},                                           &IntrinsicGeometry::computeInteriorVertexIndices);
  addQuantity(faceIndicesQ,               {},                                           &IntrinsicGeometry::computeFaceIndices);
  addQuantity(edgeIndicesQ,               {},                                           &IntrinsicGeometry::computeEdgeIndices);
  addQuantity(halfedgeIndicesQ,           {},                                           &IntrinsicGeometry::computeHalfedgeIndices);

  // == Operators
  addQuantity(basicDECOperatorsQ,         {&vertexDualAreasQ, &edgeCotanWeightsQ, &faceAreasQ, &vertexIndicesQ, &faceIndicesQ, &edgeIndicesQ},      &IntrinsicGeometry::computeBasicDECOperators);
  addQuantity(zeroFormWeakLaplacianQ,     {&basicDECOperatorsQ},                        &IntrinsicGeometry::computeZeroFormWeakLaplacian);
  // clang-format on
}

void IntrinsicGeometry::recomputeQuantities() {
  for (DependentQuantity* q : allQuantities) {
    q->computed = false;
  }
  for (DependentQuantity* q : allQuantities) {
    q->ensureHaveIfRequired();
  }
}

void IntrinsicGeometry::verifyTriangular(HalfedgeMesh* m) {
  if (!m->isSimplicial()) {
    throw std::logic_error("Only implemented for triangular meshes");
  }
}

// === Quantity implementations

 void IntrinsicGeometry::computeFaceAreas() {
  verifyTriangular(mesh);

  // TODO try these for better accuracy in near-degenerate triangles?
  // "Miscalculating Area and Angles of a Needle-like Triangle" https://www.cs.unc.edu/~snoeyink/c/c205/Triangle.pdf

  faceAreas = FaceData<double>(mesh);
  for (FacePtr f : mesh->faces()) {

    // Herons formulat
    double a = edgeLengths[f.halfedge().edge()];
    double b = edgeLengths[f.halfedge().next().edge()];
    double c = edgeLengths[f.halfedge().next().next().edge()];

    // Heron
    double s = (a + b + c) / 2.0;
    double area = std::sqrt(s * (s - a) * (s - b) * (s - c));

    faceAreas[f] = area;
  }
}


 void IntrinsicGeometry::computeVertexDualAreas() {
  vertexDualAreas = VertexData<double>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    double A = 0;
    for (FacePtr f : v.adjacentFaces()) {
      A += faceAreas[f];
    }
    vertexDualAreas[v] = A / 3.0;
  }
}

 void IntrinsicGeometry::computeHalfedgeFaceCoords() {
  verifyTriangular(mesh);
  halfedgeFaceCoords = HalfedgeData<Complex>(mesh);

  for (FacePtr f : mesh->faces()) {

    HalfedgePtr he0 = f.halfedge();
    HalfedgePtr he1 = he0.next();
    HalfedgePtr he2 = he1.next();

    // Angles measured against he0
    halfedgeFaceCoords[he0] = Complex(edgeLengths[he0.edge()], 0.0);

    // Second halfedge
    double theta1 = PI - halfedgeOppositeAngles[he2];
    halfedgeFaceCoords[he1] = std::exp(IM_I * theta1) * Complex(edgeLengths[he1.edge()], 0.0);

    // Third halfedge
    double theta2 = halfedgeOppositeAngles[he1];
    halfedgeFaceCoords[he2] = -std::exp(IM_I * theta2) * Complex(edgeLengths[he2.edge()], 0.0);
  }

  for (HalfedgePtr he : mesh->imaginaryHalfedges()) {
    halfedgeFaceCoords[he] =
        std::numeric_limits<double>::quiet_NaN(); // using this basis is never a good idea, so NaN-out
  }
}


 void IntrinsicGeometry::computeFaceTransportCoefs() {

  faceTransportCoefs = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->realHalfedges()) {
    if (he.twin().isReal()) {
      Complex angleInSource = halfedgeFaceCoords[he];
      Complex desiredAngleInTarget = -halfedgeFaceCoords[he.twin()];
      faceTransportCoefs[he] = desiredAngleInTarget / angleInSource;
    }
  }
}

 void IntrinsicGeometry::computeVertexTransportCoefs() {

  vertexTransportCoefs = HalfedgeData<Complex>(mesh);

  for (HalfedgePtr he : mesh->allHalfedges()) {
    Complex angleInSource = halfedgeVertexCoords[he];
    Complex desiredAngleInTarget = -halfedgeVertexCoords[he.twin()];
    vertexTransportCoefs[he] = desiredAngleInTarget / angleInSource;
  }
}


 void IntrinsicGeometry::computeHalfedgeOppositeAngles() {
  verifyTriangular(mesh);

  halfedgeOppositeAngles = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->allHalfedges()) {
    if (he.isReal()) {

      double lOpp = edgeLengths[he.edge()];
      double lA = edgeLengths[he.next().edge()];
      double lB = edgeLengths[he.next().next().edge()];

      double q = (lA * lA + lB * lB - lOpp * lOpp) / (2. * lA * lB);
      q = clamp(q, -1.0, 1.0);
      double angle = std::acos(q);

      halfedgeOppositeAngles[he] = angle;
    } else {
      halfedgeOppositeAngles[he] = std::numeric_limits<double>::quiet_NaN();
    }
  }
}


 void IntrinsicGeometry::computeHalfedgeCotanWeights() {
  halfedgeCotanWeights = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->realHalfedges()) {
    // halfedgeCotanWeights[he] = 1.0 / std::tan(halfedgeOppositeAngles[he]);
    halfedgeCotanWeights[he] = std::tan(PI / 2.0 - halfedgeOppositeAngles[he]);
  }
  for (HalfedgePtr he : mesh->imaginaryHalfedges()) {
    halfedgeCotanWeights[he] = std::numeric_limits<double>::quiet_NaN();
  }
}


 void IntrinsicGeometry::computeEdgeCotanWeights() {
  edgeCotanWeights = EdgeData<double>(mesh);
  for (EdgePtr e : mesh->edges()) {
    double weight = halfedgeCotanWeights[e.halfedge()];
    if (e.halfedge().twin().isReal()) {
      weight += halfedgeCotanWeights[e.halfedge().twin()];
    }
    edgeCotanWeights[e] = weight;
  }
}


 void IntrinsicGeometry::computeVertexAngleDefects() {
  verifyTriangular(mesh);

  vertexAngleDefects = VertexData<double>(mesh);
  for (VertexPtr v : mesh->vertices()) {
    if (v.isBoundary()) {
      // Convention: no curvature at boundary
      vertexAngleDefects[v] = 0;
    } else {
      double sum = 0;
      for (HalfedgePtr he : v.outgoingHalfedges()) {
        sum += halfedgeOppositeAngles[he.next()];
      }
      vertexAngleDefects[v] = 2. * PI - sum;
    }
  }
}


 void IntrinsicGeometry::computeHalfedgeRescaledOppositeAngles() {
  halfedgeRescaledOppositeAngles = HalfedgeData<double>(mesh);
  for (HalfedgePtr he : mesh->realHalfedges()) {
    double origSum = 2. * PI - vertexAngleDefects[he.next().next().vertex()];
    halfedgeRescaledOppositeAngles[he] = halfedgeOppositeAngles[he] * 2. * PI / origSum;
  }
}


 void IntrinsicGeometry::computeHalfedgeVertexCoords() {
  verifyTriangular(mesh);

  halfedgeVertexCoords = HalfedgeData<Complex>(mesh);

  for (VertexPtr v : mesh->vertices()) {

    if (v.isBoundary()) {

      // First, check what angle we associated with the boundary wedge
      // (recall that in a manifold triangle mesh, there can be at most one boundary wedge)
      double angleSum = 0;
      HalfedgePtr afterBoundaryHe;
      for (HalfedgePtr he : v.outgoingHalfedges()) {
        if (he.isReal()) {
          angleSum += halfedgeRescaledOppositeAngles[he.next()];
        }
        if (!he.twin().isReal()) {
          afterBoundaryHe = he;
        }
      }
      double boundaryAngle = 2 * PI - angleSum;

      // Now, loop like in the usual case, but substitute the boundary value when needed
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      HalfedgePtr firstHe = v.halfedge();
      HalfedgePtr currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        if (currHe.isReal()) {
          coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
          currHe = currHe.next().next().twin();
        } else {
          coordSum += boundaryAngle;
          currHe = afterBoundaryHe;
        }
      } while (currHe != firstHe);


    } else {
      double coordSum = 0.0;

      // Custom loop to orbit CCW
      HalfedgePtr firstHe = v.halfedge();
      HalfedgePtr currHe = firstHe;
      do {
        halfedgeVertexCoords[currHe] = std::exp(coordSum * IM_I);
        coordSum += halfedgeRescaledOppositeAngles[currHe.next()];
        currHe = currHe.next().next().twin();
      } while (currHe != firstHe);
    }
  }
}

 void IntrinsicGeometry::computeVertexIndices() { vertexIndices = mesh->getVertexIndices(); }

 void IntrinsicGeometry::computeInteriorVertexIndices() {
  interiorVertexIndices = mesh->getInteriorVertexIndices();
}

 void IntrinsicGeometry::computeFaceIndices() { faceIndices = mesh->getFaceIndices(); }

 void IntrinsicGeometry::computeEdgeIndices() { edgeIndices = mesh->getEdgeIndices(); }


 void IntrinsicGeometry::computeHalfedgeIndices() { halfedgeIndices = mesh->getHalfedgeIndices(); }


 void IntrinsicGeometry::computeBasicDECOperators() {

  { // Hodge 0
    size_t nVerts = mesh->nVertices();
    Eigen::VectorXd hodge0V(nVerts);
    for (VertexPtr v : mesh->vertices()) {
      double primalArea = 1.0;
      double dualArea = vertexDualAreas[v];
      double ratio = dualArea / primalArea;
      size_t iV = vertexIndices[v];
      hodge0V[iV] = ratio;
    }

    hodge0 = hodge0V.asDiagonal();
    hodge0Inv = hodge0V.asDiagonal().inverse();
  }


  { // Hodge 1
    size_t nEdges = mesh->nEdges();
    Eigen::VectorXd hodge1V(nEdges);
    for (EdgePtr e : mesh->edges()) {
      double ratio = edgeCotanWeights[e];
      size_t iE = edgeIndices[e];
      hodge1V[iE] = ratio;
    }

    hodge1 = hodge1V.asDiagonal();
    hodge1Inv = hodge1V.asDiagonal().inverse();
  }

  { // Hodge 2
    size_t nFaces = mesh->nFaces();
    Eigen::VectorXd hodge2V(nFaces);
    for (FacePtr f : mesh->faces()) {
      double primalArea = faceAreas[f];
      double dualArea = 1.0;
      double ratio = dualArea / primalArea;

      size_t iF = faceIndices[f];
      hodge2V[iF] = ratio;
    }
    hodge2 = hodge2V.asDiagonal();
    hodge2Inv = hodge2V.asDiagonal().inverse();
  }


  d0 = buildDerivative0(mesh);
  d1 = buildDerivative1(mesh);
}


 void IntrinsicGeometry::computeZeroFormWeakLaplacian() {
  zeroFormWeakLaplacian = d0.transpose() * hodge1 * d0;
}


} // namespace geometrycentral
