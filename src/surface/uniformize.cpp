#include "geometrycentral/surface/uniformize.h"


#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"


namespace geometrycentral {
namespace surface {

EdgeData<double> uniformizeDisk(IntrinsicGeometryInterface& geometry, bool withEdgeFlips) {
  HalfedgeMesh& mesh = geometry.mesh;


  // Sanity checks
  if (withEdgeFlips) {
    throw std::runtime_error("edge flips not fully implemented/tested");
  }

  // Check that it's a disk
  if (mesh.nBoundaryLoops() != 1) {
    throw std::runtime_error("parameterizeDisk(): input mesh must have exactly one boundary loop");
  }
  if (mesh.eulerCharacteristic() != 2) {
    throw std::runtime_error("parameterizeDisk(): input mesh must be topological disk");
  }


  geometry.requireEdgeLengths();

  // A geometry for the new edge lengths, which we will update
  EdgeLengthGeometry lengthGeom(mesh, geometry.edgeLengths);
  lengthGeom.requireCotanLaplacian();
  lengthGeom.requireVertexGaussianCurvatures();

  VertexData<double> u(mesh, 0.);


  // Used to decompose interior and boundary components
  size_t N = mesh.nVertices();
  size_t N_b = 0;
  Vector<bool> isInterior(N);
  for (size_t i = 0; i < N; i++) {
    isInterior(i) = !mesh.vertex(i).isBoundary();
    if (!isInterior(i)) {
      N_b++;
    }
  }


  // The four sub-blocks of the matrix are now in
  // decomp.AA, decomp.AB, decomp.BA, decomp.BB

  int nMaxIters = 10;
  for (int iIter = 0; iIter < nMaxIters; iIter++) {

    if (withEdgeFlips) {
      // TODO
    }

    // Boundary conditions (Dirichlet)
    VertexData<double> vertRHS(mesh, 0.);
    for (Vertex v : mesh.vertices()) {
      if (!v.isBoundary()) {
        vertRHS[v] = -lengthGeom.vertexGaussianCurvatures[v];
      }
    }

    // Build the cot-Laplace operator
    Vector<double> rhsVals = vertRHS.toVector();
    Vector<double> bcVals = Vector<double>::Zero(N_b); // "natural" dirichlet boundary conditions

    // Step the flow
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(lengthGeom.cotanLaplacian, isInterior, true);

    // Split up the rhs vector
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, rhsVals, rhsValsA, rhsValsB);

    // Solve problem
    Vector<double> combinedRHS = rhsValsA - decomp.AB * bcVals;
    Vector<double> Aresult = solve(decomp.AA, combinedRHS);

    // Combine the two boundary conditions and interior solution to a full vector
    Vector<double> result = reassembleVector(decomp, Aresult, bcVals);

    // Update edge lengths
    u.fromVector(result);
    double maxRelChange = 0;
    for (Edge e : mesh.edges()) {
      double u1 = u[e.halfedge().vertex()];
      double u2 = u[e.halfedge().twin().vertex()];
      double s = std::exp((u1 + u2) / 2.);
      double oldLen = lengthGeom.inputEdgeLengths[e];
      double newLen = oldLen * s;

      //std::cout << "u1 = " << u1 << " u2 = " << u2 << " s = " << s << " oldLen = " << oldLen << " newLen = " << newLen
                //<< std::endl;

      double relChange = std::fabs(newLen - oldLen) / oldLen;
      maxRelChange = std::fmax(relChange, maxRelChange);

      lengthGeom.inputEdgeLengths[e] = newLen;
    }

    // TODO convergence check?
    std::cout << " uniformize max change = " << maxRelChange << std::endl;

    lengthGeom.refreshQuantities();
  }

  geometry.unrequireEdgeLengths();

  return lengthGeom.inputEdgeLengths;
}

} // namespace surface
} // namespace geometrycentral
