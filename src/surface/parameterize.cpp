#include "geometrycentral/surface/parameterize.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/uniformize.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace geometrycentral {
namespace surface {

VertexData<Vector2> parameterizeDisk(ManifoldSurfaceMesh& origMesh, IntrinsicGeometryInterface& origGeom) {
  // Check that it's a (punctured) disk

  /*
   if ((long long int)origMesh.eulerCharacteristic() - (long long int)(2 * origMesh.nBoundaryLoops()) != -2) {
    long long int val =
        ((long long int)origMesh.eulerCharacteristic() - (long long int)(2 * origMesh.nBoundaryLoops()));
    throw std::runtime_error("parameterizeDisk(): input origMesh must be a (possibly punctured) disk, chi - 2b = " +
                             std::to_string(val));
  }
  */

  // Get uniformized edge lengths
  // Copy the mesh, since we will flip its edges
  std::unique_ptr<ManifoldSurfaceMesh> meshPtr = origMesh.copy();
  ManifoldSurfaceMesh& mesh = *meshPtr;
  origGeom.requireEdgeLengths();
  EdgeData<double> copyLens = origGeom.edgeLengths.reinterpretTo(mesh);
  EdgeLengthGeometry geometry(mesh, copyLens);
  origGeom.unrequireEdgeLengths();
  EdgeData<double> uLens = uniformizeDisk(mesh, geometry, true);

  // Layout
  VertexData<Vector2> coords(mesh);
  VertexData<char> haveCoords(mesh, false);

  // == Layout until finished

  // <n_neighbors, -n_round, vert>
  // prefer verts with more neighbors, breaking ties by earlier verts
  using WeightedVertex = std::tuple<int, int, Vertex>;
  std::priority_queue<WeightedVertex> pq;
  int layoutRound = 0;

  // Helper to add vertices to the layout queue if they have enough neighbors to layout
  auto considerVertex = [&](Vertex v) {
    if (haveCoords[v]) return;

    int neighCount = 0;
    for (Halfedge he : v.outgoingHalfedges()) {
      if (!he.isInterior()) continue;
      Halfedge heOpp = he.next();
      if (haveCoords[heOpp.vertex()] && haveCoords[heOpp.twin().vertex()]) neighCount++;
    }

    if (neighCount >= 1) pq.emplace(neighCount, -layoutRound, v);
  };

  // initialize the layout
  Halfedge he0 = mesh.halfedge(0);
  Vertex v0 = he0.vertex();
  Vertex v1 = he0.twin().vertex();
  coords[v0] = Vector2{0., 0.};
  coords[v1] = Vector2{uLens[he0.edge()], 0.};
  haveCoords[v0] = true;
  haveCoords[v1] = true;
  considerVertex(he0.next().next().vertex());
  considerVertex(he0.twin().next().next().vertex());

  // layout until finished
  while (!pq.empty()) {

    int topCount = std::get<0>(pq.top());
    Vertex topVert = std::get<2>(pq.top());
    pq.pop();

    if (haveCoords[topVert]) continue;

    // std::cout << "laying out " << topVert << " with " << topCount << " neighbors" << std::endl;

    // Compute a new position for the vertex, as an average of laid out position from all neighbors
    Vector2 avgPos{0., 0.};
    int avgCount = 0;
    for (Halfedge he : topVert.outgoingHalfedges()) {
      if (!he.isInterior()) continue;
      Halfedge heOpp = he.next();
      Vertex vA = heOpp.vertex();
      Vertex vB = heOpp.twin().vertex();
      if (haveCoords[vA] && haveCoords[vB]) {

        Vector2 pA = coords[vA];
        Vector2 pB = coords[vB];
        double lCA = uLens[he.edge()];
        double lBC = uLens[heOpp.next().edge()];

        Vector2 newPos = layoutTriangleVertex(pA, pB, lBC, lCA);
        avgPos += newPos;
        avgCount++;
      }
    }

    // set the new positions
    coords[topVert] = avgPos / avgCount;
    haveCoords[topVert] = true;

    // add neighbors
    for (Vertex vn : topVert.adjacentVertices()) {
      considerVertex(vn);
    }

    layoutRound++;
  }

  return coords.reinterpretTo(origMesh);
}

VertexData<Vector3> parameterizeSphere(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry) {

  // Check that it's a sphere
  if (mesh.eulerCharacteristic() != -2 ||
      mesh.nBoundaryLoops() != 0) {
     throw std::runtime_error("parameterizeSphere(): mesh must have spherical topology (and no boundary)");
  }
  
  // Check that it's triangulated
  for( Face f : mesh.faces() ) {
     if( f.degree() != 3 ) {
        throw std::runtime_error("parameterizeSphere(): all mesh faces must be triangles");
     }
  }

  VertexData<Vector3> coords(mesh);

  // Select a triangle to pin, and compute locations for its three vertices
  Face fp = mesh.face(0);
  array<Vertex,3> vp; // pinned vertices
  VertexData<Vector2> param( mesh ); // 2D parameterization
  Halfedge h = fp.halfedge();
  for( int k = 0; k < 3; k++ ) {
     // grab vertices of the pinned triangle
     vp[k] = h.vertex();
     h = h.next();

     // construct vertices of an equilateral triangle
     double theta = 2.*k*M_PI/3.;
     param[vp[k]] = Vector2{ cos(theta), sin(theta) };
  }

  // Assign indices to non-pinned vertices, and flag pinned vertices
  VertexData<size_t> index;
  size_t nV = 0;
  size_t pinned = -1;
  for( Vertex v : mesh.vertices() ) {
     if( v == vp[0] && v == vp[1] && v == vp[2] ) {
        index[v] = pinned;
     } else {
        index[v] = nV++;
     }
  }

  // Build (weak) cotan-Laplace operator with three vertices pinned
  int n = mesh.nVertices() - 3;
  SparseMatrix<double> L( n, n );
  std::vector<Vector2> b( n, Vector2{0.,0.} );
  typedef Eigen::triplet<double> Entry;
  std::vector<Entry> entries;

  // set entries by iterating over halfedges
  geometry.requireHalfedgeCotanWeights();
  for( Halfedge h : mesh.halfedges() ) {
     // get the endpoints and their indices
     Vertex vi = h.vertex();
     Vertex vj = h.vertex().twin();
     size_t i = index[vi];
     size_t j = index[vj];

     // get the (half)edge weight
     double cotTheta = halfedgeCotanWeights[h];

     if( i != pinned ) {
        entries.push_back( Triplet( i, i, cotTheta ));
        if( j != pinned ) {
           entries.push_back( Triplet( i, j, -cotTheta ));
        }
        else {
           b(i) += cotTheta * param[vj];
        }
     }

     if( j != pinned ) {
        entries.push_back( Triplet( j, j, cotTheta ));
        if( i != pinned ) {
           entries.push_back( Triplet( j, i, -cotTheta ));
        }
        else {
           b(j) += cotTheta * param[vi];
        }
     }
  }

  // solve for the two coordinate functions
  PositiveDefiniteSolver<double> solver( L );
  for( int k = 0; k < 2; j++ ) {

     // build a column vector for the kth component of the right-hand side
     Vector bk( n );
     for( size_t i = 0; i < n; i++ ) {
        bk[i] = b[i][k];
     }

     // solve for the kth coordinate function in the plane
     Vector<double> xk = solver.solve( bk );

     // copy into 2D coordinate vectors
     for( Vertex vi : mesh.vertices() ) {
        size_t i = index[vi];
        if( i != pinned ) {
           param[vi][k] = xk[i];
        }
     }
  }

  // stereographically map from equilateral triangle to the sphere
  const double scaleFactor = 1e5;
  for( Vertex v : mesh.vertices() ) {
     Vector3 x&( coords[v] );
     Vector2 X = param[v];

     // make the boundary equilateral triangle very big, so
     // that its complement (corresponding to the removed
     // triangle projects to a very small area on the unit sphere)
     X *= scaleFactor;

     // apply stereographic projection
     x = Vector3 {
         2.*X[0],
         2.*X[1],
         X[0]*X[0] + X[1]*X[1] - 1.
     } / X[0]*X[0] + X[1]*X[1] + 1.;
  }

  // TODO Möbius balancing
}

} // namespace surface
} // namespace geometrycentral
