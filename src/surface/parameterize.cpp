#include "geometrycentral/surface/parameterize.h"

#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/uniformize.h"
#include "geometrycentral/utilities/elementary_geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

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

Eigen::Vector3d toEigen( const Vector3& u ) {
   Eigen::Vector3d v;
   v << u[0], u[1], u[2];
   return v;
}

Vector3 toVector3( const Eigen::Vector3d& u ) {
   return Vector3{ u(0), u(1), u(2) };
}

VertexData<Vector3> parameterizeSphere(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geometry) {

   // Compute a conformal parameterization to the sphere, in three steps:
   //
   //    1. Map the surface minus a triangle into an equilateral triangle in the plane.
   //    2. Stereographically project from the plane to the sphere.
   //    3. Compute the unique Möbius transformation that puts the center of mass at the origin.
   //
   // This algorithm basically builds on two papers:
   //
   //   A Linear Variational Principle for Riemann Mappings and Discrete Conformality
   //   Nadav Dym, Raz Slutsky, Yaron Lipman
   //   Proceedings of the National Academy of Sciences (PNAS), 116(3), 2019
   //
   //   Möbius Registration
   //   Alex Baden, Keenan Crane, Misha Kazhdan
   //   Symposium on Geometry Processing (SGP), 37(5), 2018
   //
   // The first paper basically says that if you map a surface with boundary into an equilateral
   // triangle via a harmonic map, but let the boundary vertices "slide" along the three triangle
   // edges (which is just a linear constraint), then you'll get something that converges to a
   // conformal map, rather than just a harmonic one.  This fact is essentially an application of
   // Eells & Wood 1975, which likewise says that harmonic maps between surfaces with appropriate
   // Euler characteristic are automatically conformal.  The use of an equilateral triangle can
   // also be justified from an orbifold perspective.  In this algorithm, we consider a very
   // special case where the boundary is just a single triangle.  Hence, to compute the harmonic
   // map we just need to compute a discrete harmonic map with Dirichlet conditions at three
   // vertices, using the usual cotan matrix L:
   //
   //    Lx = 0 on V \ {i,j,k}
   //    x_i = (1,0)
   //    x_j = (-1/2, √3/2)
   //    x_k = (-1/2,-√3/2)
   //
   // The second paper accounts for the fact that an arbitrary stereographic projection can
   // result in extreme distortion; it gives a simple procedure for computing the unique
   // Möbius transformation that puts the center of mass (relative to the density on the original
   // surface) at the origin, giving the "best possible" map from the 2D domain to the sphere.
   // See Algorithm 1 of Baden et al for pseudocode.  The code below makes a couple refinements,
   // not described in the original paper, which makes this process faster and more numerically
   // stable, namely:
   //
   //    - It picks an initial scaling for the 2D domain so that the removed triangle has
   //      approximately the same area (proportionally) as in the original mesh, helping
   //      to provide a good initialization.
   //    - It also uses a more intelligent step size, essentially by viewing the differential
   //      of the energy as a tangent vector to the space of Möbius transformations, and
   //      using the exponential map (relative to the metric in the Poincaré ball model) to
   //      compute the actual update.
   //
   // As in the original method, this procedure is guaranteed to converge to the optimal
   // solution, and typically does so using fewer than 10 iterations on a small optimization
   // problem with only three degrees of freedom.
   //

   // TODO Should probably check that we have a pure triangle mesh of a sphere.

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
   }

   // TODO refactor into separate subroutine
   // Compute a canonical Möbius transformation, a la
   // Baden et al, "Möbius Registration" (SGP 2018)
   const double eps = 1e-7; // stopping tolerance
   const int maxIter = 100; // if it takes more than a hundred iterations, something is wrong anyway!
   VertexPositionGeometry sphereGeometry( mesh, coords );
   for( int iter = 0; iter < maxIter; iter++ )
   {
      // compute center of mass µ and Jacobian J
      Vector3 mu{ 0., 0., 0. };
      Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
      Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
      for( Face tau : mesh.faces() )
      {
         double Atau = geometry.faceAreas[tau];
         Vector3 Ctau = (sphereGeometry.faceCentroid(tau)).unit();
         Eigen::Vector3d C = toEigen(Ctau);

         mu += Atau * Ctau;
         J += 2.*Atau * ( I - C*C.transpose() );
      }

      // check for convergence
      std::cerr << "|µ|: " << mu.norm() << std::endl;
      if( mu.norm() < eps ) {
         break;
      }

      // compute inversion center c = -J⁻¹μ
      Vector3 c = -toVector3( J.inverse() * toEigen(mu) );
      c = c * tanh( std::min( c.norm(), 2. ))/c.norm(); // stupid hyperbolic tricks...

      // apply inversion
      for( Vertex i : mesh.vertices() ) {
         Vector3& v = sphereGeometry.vertexPositions[i];
         v = (1-c.norm2()) * (v+c)/(v+c).norm2() + c;
      }
   }

   return sphereGeometry.vertexPositions;
}

} // namespace surface
} // namespace geometrycentral
