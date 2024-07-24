#pragma once

#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {

struct EmbedConvexResult {
  std::unique_ptr<ManifoldSurfaceMesh> mesh;
  std::unique_ptr<VertexPositionGeometry> geometry;
  bool success = false; // obtained an isometric embedding
};

struct EmbedConvexOptions {
  double initialStepSize = 1.;            // initial step size for optimization
  int maxSteps = 20;                      // maximum number of optimization steps
  int maxLineSearchSteps = 16;            // maximum number of line search steps
  double embeddingTolerance = 1e-3;       // maximum angle defect for any radial edge
  int maxNewtonIterations = 20;           // maximum number of steps for inner Newton solver
  int maxNewtonLineSearchSteps = 32;      // maximum number of line search steps for inner Newton solver
  double newtonTolerance = 1e-4;          // l2 norm of residual for convergence of Newton's method
  double metricConvexityTolerance = 1e-5; // how negative do we allow input Gaussian curvature to be?
  double edgeConvexityTolerance = 1e-3;   // how negative do we allow mean curvature to be?
  bool verbose = false;                   // whether to display diagnostic output
};
extern const EmbedConvexOptions defaultEmbedConvexOptions;

// Find an embedding given lengths on edges
EmbedConvexResult embedConvex(ManifoldSurfaceMesh& triangulation, EdgeData<double>& edgeLengths,
                              EmbedConvexOptions options = defaultEmbedConvexOptions);

// Find an embedding given 2D texture coordinates that define the
// intrinsic edge lengths.  NOTE: lengths that are not compatible
// across seams in the UV map will be averaged.
EmbedConvexResult embedConvex(ManifoldSurfaceMesh& triangulation, CornerData<Vector2>& textureCoordinates,
                              EmbedConvexOptions options = defaultEmbedConvexOptions);

// Find an embedding given an intrinsic triangulation
EmbedConvexResult embedConvex(ManifoldSurfaceMesh& triangulation, IntrinsicGeometryInterface& intrinsicGeometry,
                              EmbedConvexOptions options = defaultEmbedConvexOptions);

class ConvexEmbedder {
public:
  ConvexEmbedder(ManifoldSurfaceMesh& trianguation, EdgeData<double>& edgeLengths,
                 EmbedConvexOptions options = defaultEmbedConvexOptions);
  ConvexEmbedder(ManifoldSurfaceMesh& trianguation, CornerData<Vector2>& textureCoordinates,
                 EmbedConvexOptions options = defaultEmbedConvexOptions);
  ConvexEmbedder(ManifoldSurfaceMesh& trianguation, IntrinsicGeometryInterface& intrinsicGeometry,
                 EmbedConvexOptions options = defaultEmbedConvexOptions);
  ~ConvexEmbedder();

  // Algorithm parameters
  EmbedConvexOptions options;

  // Initialize the embedding, returning false if the
  // given triangulation and metric do not meet the
  // criteria for convex embedding (no boundary,
  // triangulated, genus zero, positive angle defect).
  bool init();

  // Embed the given metric, returning true if embedding was successful.
  // Note that this method "fails gracefully" in the sense that if we
  // don't find an embedding that meets the specified tolerance, we still
  // construct vertex positions that are as close as we could get, by
  // averaging the most recent corner coordinates.
  bool embed();

  // Take a single embedding step (for visualization/debugging)
  void stepEmbedding();

  // Compute extrinsic vertex positions at corners from the
  // current intrinsic embedding, stored in `localLayout`.
  // Also compute average coordinates at vertices, stored
  // in `embeding`.
  void refreshVertexCoordinates();

  // Original metric
  ManifoldSurfaceMesh& originalMesh;
  IntrinsicGeometryInterface* originalMetric;

  // Intrinsic triangulation used for the embedding
  IntegerCoordinatesIntrinsicTriangulation* intrinsicTriangulation;
  ManifoldSurfaceMesh* currentMesh; // points to intrinsicTriangulation->intrinsicMesh

  // Embedding
  VertexData<double> r;            // length of radial edges connecting each vertex to the central apex
  CornerData<Vector3> localLayout; // coordinates in R^3 for each corner of currentMesh (may be discontinuous)
  VertexData<Vector3> embedding; // coordinates in R^3 for each vertex of currentMesh (average of discontinuous values)
  Vector3 apex;                  // location of central vertex, making tetrahedra with each triangle

protected:
  bool isMetricConvex();    // do the input edge lengths exhibit vertex angle sums ≤ 2π?
  bool isEmbeddingConvex(); // does any boundary edge have negative curvature?

  double hessianOffDiagonal(Halfedge ij);
  double hessianDiagonal(Vertex i);
  std::pair<SparseMatrix<double>, bool> buildHessian();

  // Returns the minimum residual of the triangle
  // inequality at any corner; if this value is
  // negative, the triangle inequality has been violated.
  double minTriSlack();

  // Local geometry
  double l(Edge ij);        // length of boundary edge ij
  double l(Halfedge ij);    // length of boundary edge ij.edge()
  double alpha(Halfedge h); // dihedral angles around boundary edges
  double beta(Corner c);    // dihedral angles around radial edges
  double gamma(Corner ijk); // interior angle of boundary triangle ijk at corner i
  double theta(Edge ij);    // sum of dihedral angles around boundary edges
  double omega(Vertex i);   // sum of dihedral angles around interior edges
  double delta(Vertex i);   // angle defects at boundary vertices
  double kappa(Vertex i);   // curvature around interior edges
  double phi(Halfedge ij);  // interior angle at vertex a of triangle aij (where a is the apex)
  double rho(Halfedge ij);  // angle between boundary edge ij and the radius aj
  double circumradius(Face f);

  // These routines directly compute the cos/sin
  // of an angle, without trigonometric functions.
  double cosRho(Halfedge ij);                       // returns cosine of rho(ij)
  double cosGamma(Corner ij);                       // returns cosine of gamma(ijk)
  std::pair<double, double> cosSinRho(Halfedge ij); // returns cosine and sine of rho(ij)
  std::pair<double, double> cosSinGamma(Corner ij); // returns cosine and sine of gamma(ijk)

  // Utility functions
  double cot(double angle);
  double acosClamp(double angle);
  double interiorAngleCosine(double a, double b, double c);
  double interiorAngle(double a, double b, double c);

  // Initialization
  double initialRadius();
  bool initializeRadii();

  // Optimization
  bool solveNewton(VertexData<double> kappaStar);
  bool embeddingConverged(); // did optimization converge to the given tolerance?

  // Edge flipping
  void flipToConvex(); // try to flip all edges with negative mean curvature

  // Keep track of whether originalMetric was locally
  // allocated, or points to an external instance
  bool originalMetricLocallyAllocated;
};

} // namespace surface
} // namespace geometrycentral
