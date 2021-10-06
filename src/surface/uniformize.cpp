#include "geometrycentral/surface/uniformize.h"


#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/simple_idt.h"


namespace geometrycentral {
namespace surface {

EdgeData<double> uniformizeDisk(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry, bool withEdgeFlips) {

  // Check that it's a (punctured) disk
  /*
  if ((long long int)mesh.eulerCharacteristic() - (long long int)(2 * mesh.nBoundaryLoops()) != -2) {
    long long int val = ((long long int)mesh.eulerCharacteristic() - (long long int)(2 * mesh.nBoundaryLoops()));
    throw std::runtime_error("uniformizeDisk(): input must be a (possibly punctured) disk, chi - 2b = " +
                             std::to_string(val));
  }
  */


  geometry.requireEdgeLengths();

  // A geometry for the new edge lengths, which we will update
  EdgeLengthGeometry lengthGeom(mesh, geometry.edgeLengths);
  if (withEdgeFlips) {
    flipToDelaunay(mesh, lengthGeom.edgeLengths, FlipType::Hyperbolic);
  }
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


  int nMaxIters = 50;
  for (int iIter = 0; iIter < nMaxIters; iIter++) {

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
    SparseMatrix<double> A = lengthGeom.cotanLaplacian;
    BlockDecompositionResult<double> decomp = blockDecomposeSquare(A, isInterior, true);

    // Split up the rhs vector
    Vector<double> rhsValsA, rhsValsB;
    decomposeVector(decomp, rhsVals, rhsValsA, rhsValsB);

    // Solve problem
    Vector<double> combinedRHS = rhsValsA - decomp.AB * bcVals;
    Vector<double> Aresult = 0.5 * solve(decomp.AA, combinedRHS);

    // Combine the two boundary conditions and interior solution to a full vector
    Vector<double> result = reassembleVector(decomp, Aresult, bcVals);

    // Update edge lengths
    u.fromVector(result);
    double maxRelChange = 0;
    for (Edge e : mesh.edges()) {
      double u1 = u[e.halfedge().vertex()];
      double u2 = u[e.halfedge().twin().vertex()];
      double s = std::exp((u1 + u2) / 2.);
      double oldLen = lengthGeom.edgeLengths[e];
      double newLen = oldLen * s;

      // std::cout << "u1 = " << u1 << " u2 = " << u2 << " s = " << s << " oldLen = " << oldLen << " newLen = " <<
      // newLen
      //<< std::endl;

      double relChange = std::fabs(newLen - oldLen) / oldLen;
      maxRelChange = std::fmax(relChange, maxRelChange);

      lengthGeom.edgeLengths[e] = newLen;
    }

    // TODO convergence check?
    std::cout << " uniformize max change = " << maxRelChange << std::endl;

    if (withEdgeFlips) {
      flipToDelaunay(mesh, lengthGeom.edgeLengths, FlipType::Hyperbolic);
    }
    lengthGeom.refreshQuantities();
  }

  geometry.unrequireEdgeLengths();

  return lengthGeom.edgeLengths;
}

} // namespace surface
} // namespace geometrycentral
