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

  // // Check that it's a sphere
  // if (mesh.eulerCharacteristic() != -2 ||
  //     mesh.nBoundaryLoops() != 0) {
  //    throw std::runtime_error("parameterizeSphere(): mesh must have spherical topology (and no boundary)");
  // }
  // 
  // // Check that it's triangulated
  // for( Face f : mesh.faces() ) {
  //    if( f.degree() != 3 ) {
  //       throw std::runtime_error("parameterizeSphere(): all mesh faces must be triangles");
  //    }
  // }

  VertexData<Vector3> coords(mesh);

  // Select a triangle to pin, and compute locations for its three vertices
  Face fp = mesh.face(0);
  std::array<Vertex,3> vp; // pinned vertices
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
  VertexData<size_t> index( mesh );
  size_t nV = 0;
  size_t pinned = -1;
  for( Vertex v : mesh.vertices() ) {
     if( v == vp[0] || v == vp[1] || v == vp[2] ) {
        index[v] = pinned;
     } else {
        index[v] = nV;
        nV++;
     }
  }

  // Build (weak) cotan-Laplace operator with three vertices pinned
  size_t n = mesh.nVertices() - 3;
  std::vector<Vector2> b( n, Vector2{0.,0.} );
  typedef Eigen::Triplet<double> Entry;
  std::vector<Entry> entries;

  // set entries by iterating over halfedges
  geometry.requireHalfedgeCotanWeights();
  for( Halfedge h : mesh.halfedges() ) {
     // get the endpoints and their indices
     Vertex vi = h.vertex();
     Vertex vj = h.twin().vertex();
     size_t i = index[vi];
     size_t j = index[vj];

     // get the (half)edge weight
     double cotTheta = geometry.halfedgeCotanWeights[h];

     // set the entries
     if( i != pinned ) {
        entries.emplace_back( Entry( i, i, cotTheta ));
        // std::cerr << "{" << 1+i << "," << 1+i << "} -> " << cotTheta << ", ";
        if( j != pinned ) {
           entries.emplace_back( Entry( i, j, -cotTheta ));
           // std::cerr << "{" << 1+i << "," << 1+j << "} -> " << -cotTheta << ", ";
        }
        else {
           b[i] += cotTheta * param[vj];
        }
     }

     if( j != pinned ) {
        entries.emplace_back( Entry( j, j, cotTheta ));
        // std::cerr << "{" << 1+j << "," << 1+j << "} -> " << cotTheta << ", ";
        if( i != pinned ) {
           entries.emplace_back( Entry( j, i, -cotTheta ));
           // std::cerr << "{" << 1+j << "," << 1+i << "} -> " << -cotTheta << ", ";
        }
        else {
           b[j] += cotTheta * param[vi];
        }
     }
  }

  // set from triplets
  std::cerr << "set from triplets" << std::endl;
  SparseMatrix<double> L( n, n );
  L.setFromTriplets( entries.begin(), entries.end() );

  // solve for the two coordinate functions
  std::cerr << "solve for the two coordinate functions" << std::endl;
  PositiveDefiniteSolver<double> solver( L );
  for( int k = 0; k < 2; k++ ) {

     // build a column vector for the kth component of the right-hand side
     Vector<double> bk( n );
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

  // Compute scale factor needed for next step, where we use
  // stereographic projection to map the triangular domain to
  // the unit sphere.  If the image of the outer triangle is
  // too small or too big, this map will be badly distorted.
  // So, as a heuristic, we'll first scale the 2D domain so
  // that the image triangle has area roughly proportional to
  // its area in the original mesh.  Specifically, suppose we
  // apply a scale factor h; then the vertex (1,0,0) gets
  // scaled up to (h,0,0), and its x-coordinate gets mapped
  // by stereographic projection to r := 2h/(h^2 + 1).  As a
  // heuristic, we will solve for h such that πr^2 = 4π A0,
  // where A0 is the area fraction and 4π is the area of the
  // unit sphere.
  geometry.requireFaceAreas();
  double surfaceArea = 0.;
  for( Face f : mesh.faces() ) {
     surfaceArea += geometry.faceAreas[f];
  }
  double A0 = geometry.faceAreas[fp] / surfaceArea;
  double h0 = (1.-sqrt(1.-4.*A0))/(2.*sqrt(A0));
  double h1 = (1.+sqrt(1.-4.*A0))/(2.*sqrt(A0));
  double scaleFactor = std::max(h0,h1);

  // stereographically map from equilateral triangle to the sphere
  for( Vertex v : mesh.vertices() ) {
     Vector2 X = param[v];
     Vector3& x( coords[v] );

     // Make the boundary equilateral triangle bigger, so
     // that its complement (corresponding to the removed
     // triangle projects to a very small area on the unit sphere).
     // As a heuristic, we'll try to roughly match the size of the
     // outer triangle so it's proportional to its original area.
     X *= scaleFactor;

     // apply stereographic projection
     x = Vector3 {
          2.*X[0],
          2.*X[1],
          X[0]*X[0] + X[1]*X[1] - 1.
     } /( X[0]*X[0] + X[1]*X[1] + 1. );

     // x = Vector3{ X[0], X[1], 0. };
  }

  // TODO Möbius balancing

  return coords;
}

} // namespace surface
} // namespace geometrycentral
