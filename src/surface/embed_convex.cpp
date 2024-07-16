#include "geometrycentral/surface/embed_convex.h"
#include "geometrycentral/surface/simple_idt.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <queue>
using namespace std;

namespace geometrycentral {
namespace surface {

// The default options
const EmbedConvexOptions defaultEmbedConvexOptions;

EmbedConvexResult embedConvex( ManifoldSurfaceMesh& triangulation,
                               EdgeData<double>& edgeLengths,
                               EmbedConvexOptions options ) {
   EmbedConvexResult result;
   ConvexEmbedder embedder( triangulation, edgeLengths, options );

   if( !embedder.init() ) {
      result.success = false;
      return result;
   }

   result.success = embedder.embed();
   embedder.refreshVertexCoordinates();
   result.mesh = embedder.imesh->copy();
   result.geometry = unique_ptr<VertexPositionGeometry>(
         new VertexPositionGeometry(
            *result.mesh,
            embedder.embedding.reinterpretTo( *result.mesh )
            )
         );

   return result;
}

EmbedConvexResult embedConvex( ManifoldSurfaceMesh& triangulation,
                               CornerData<Vector2>& uv,
                               EmbedConvexOptions options )
{
   EmbedConvexResult result;
   ConvexEmbedder embedder( triangulation, uv, options );

   if( !embedder.init() ) {
      result.success = false;
      return result;
   }

   result.success = embedder.embed();
   embedder.refreshVertexCoordinates();
   result.mesh = embedder.imesh->copy();
   result.geometry = unique_ptr<VertexPositionGeometry>(
         new VertexPositionGeometry(
            *result.mesh,
            embedder.embedding.reinterpretTo( *result.mesh )
            )
         );

   return result;
}

EmbedConvexResult embedConvex( ManifoldSurfaceMesh& triangulation,
                               IntrinsicGeometryInterface& intrinsicGeometry,
                               EmbedConvexOptions options )
{
   EmbedConvexResult result;
   ConvexEmbedder embedder( triangulation, intrinsicGeometry, options );

   if( !embedder.init() ) {
      result.success = false;
      return result;
   }

   result.success = embedder.embed();
   embedder.refreshVertexCoordinates();
   result.mesh = embedder.imesh->copy();
   result.geometry = unique_ptr<VertexPositionGeometry>(
         new VertexPositionGeometry(
            *result.mesh,
            embedder.embedding.reinterpretTo( *result.mesh )
            )
         );

   return result;
}

ConvexEmbedder::ConvexEmbedder(
      ManifoldSurfaceMesh& triangulation,
      EdgeData<double>& edgeLengths,
      EmbedConvexOptions embedderOptions )
: originalMesh( triangulation )
{
   originalMetric = new EdgeLengthGeometry( originalMesh, edgeLengths );
   originalMetricLocallyAllocated = true;

   igeom = new IntegerCoordinatesIntrinsicTriangulation( originalMesh, *originalMetric );
   imesh = igeom->intrinsicMesh.get();

   options = embedderOptions;
}

ConvexEmbedder::ConvexEmbedder(
      ManifoldSurfaceMesh& triangulation,
      CornerData<Vector2>& textureCoordinates,
      EmbedConvexOptions embedderOptions )
: originalMesh( triangulation )
{
   // measure average length of each edge in UV space
   EdgeData<double> edgeLengths( originalMesh );
   for( Edge ij : originalMesh.edges() )
   {
      Corner ijk = ij.halfedge().corner();
      Corner jki = ij.halfedge().next().corner();
      Corner jil = ij.halfedge().twin().corner();
      Corner ilj = ij.halfedge().twin().next().corner();
      double lij = ( textureCoordinates[ijk] - textureCoordinates[jki] ).norm();
      double lji = ( textureCoordinates[jil] - textureCoordinates[ilj] ).norm();
      edgeLengths[ij] = (lij+lji)/2.;
   }

   originalMetric = new EdgeLengthGeometry( originalMesh, edgeLengths );
   originalMetricLocallyAllocated = true;

   igeom = new IntegerCoordinatesIntrinsicTriangulation( originalMesh, *originalMetric );
   imesh = igeom->intrinsicMesh.get();

   options = embedderOptions;
}

ConvexEmbedder::ConvexEmbedder(
      ManifoldSurfaceMesh& triangulation,
      IntrinsicGeometryInterface& intrinsicGeometry,
      EmbedConvexOptions embedderOptions )
   : originalMesh( triangulation )
{
   originalMetric = &intrinsicGeometry;
   originalMetricLocallyAllocated = false;

   igeom = new IntegerCoordinatesIntrinsicTriangulation( originalMesh, *originalMetric );
   imesh = igeom->intrinsicMesh.get();

   options = embedderOptions;
}

ConvexEmbedder::~ConvexEmbedder()
{
   if( originalMetricLocallyAllocated )
   {
      delete originalMetric;
   }
   delete igeom;
}

bool ConvexEmbedder::init()
{
   if( imesh->hasBoundary() ) {
      cerr << "ConvexEmbedder: cannot embed surfaces with boundary." << endl;
      return false;
   }
   if( !imesh->isTriangular() ) {
      cerr << "ConvexEmbedder: all faces of input mesh must be triangular." << endl;
      return false;
   }
   if( imesh->genus() != 0 ) {
      cerr << "ConvexEmbedder: input mesh must have spherical topology." << endl;
      return false;
   }
   if( !isMetricConvex() ) {
      cerr << "ConvexEmbedder: input metric must have vertex angle sums ≤ 2π." << endl;
      return false;
   }

   flipToDelaunay(*imesh, igeom->edgeLengths);
   if( !initializeRadii() )
   {
      cerr << "ConvexEmbedder: failed to initialize embedding." << endl;
      return false;
   }

   return true;
}

bool ConvexEmbedder::embed()
{
   // Try to embed the surface
   int step;
   for( step = 0; step < options.maxSteps; step++ )
   {
      stepEmbedding();
      if( embeddingConverged() ) {
         cout << "ConvexEmbedder: converged after " << step << " steps." << endl;
         break;
      }
   }

   if( step == options.maxSteps ) {
      cerr << "ConvexEmbedder: failed to converge to tolerance of ";
      cerr << options.tolerance << " after ";
      cerr << options.maxSteps << " steps." << endl;
      cerr << "(Try using looser tolerance or more steps?)" << endl;
      return false;
   }
   
   return true; // embedding succeeded
}

bool ConvexEmbedder::embeddingConverged() {

   // Check if the dihedral angle defect around all radial
   // edges is sufficiently close to zero (i.e., close to
   // being flat, and hence embeddable in ℝ³).
   for( Vertex i : imesh->vertices() )
   {
      if( abs(kappa(i)) > options.tolerance ) {
         return false;
      }
   }
   return true;
}

// Check if the surface metric defined by the input
// edge lengths l is convex by seeing whether every
// vertex has nonnegative angle defect
bool ConvexEmbedder::isMetricConvex() {
   const double eps = 1e-5; // tolerance for geometric checks of convexity
   for( Vertex i : imesh->vertices() ) {
      if( delta(i) < -eps ) {
         return false;
      }
   }
   return true;
}

// Check if the generalized polyhedral metric defined
// by the current radii r is convex by seeing if any
// boundary edge is concave
bool ConvexEmbedder::isEmbeddingConvex() {
   const double eps = 1e-3; // numerical tolerance
   for( Edge e : imesh->edges() ) {
      if( M_PI - theta(e) < -eps ) {
         return false;
      }
   }
   return true;
}

// The smallest possible radius R is the maximum
// circumradius of any triangle, since any shorter R
// will not define a valid tetrahedron above the
// maximal-radius triangle.
double ConvexEmbedder::initialRadius()
{
   double R = 0.;
   for( Face ijk : imesh->faces() )
   {
      double Rijk = circumradius( ijk );
      R = max( R, Rijk );
   }
   return R;
}

// Find a sufficiently large initial radius R so
// that all tetrahedra are well-defined and the
// boundary (mean) curvature is positive.
bool ConvexEmbedder::initializeRadii()
{
   // initial radii for a convex metric
   r = VertexData<double>( *imesh );
   double R = initialRadius();
   const double initRadiusFactor = 1.5;
   bool polytopeIsValid = false;
   while( true ) {
      for (Vertex i : imesh->vertices()) {
         r[i] = R;
      }

      polytopeIsValid = isEmbeddingConvex();
      if( !polytopeIsValid ) {
         cerr << "Mean curvature is negative for radius " << R << endl;
      }

      for( Vertex i : imesh->vertices() ) {
         if( kappa(i) < -1e-3 || delta(i) < kappa(i) ) {
            polytopeIsValid = false;
            cerr << "Sectional curvature " << kappa(i) << " is out of bounds for radius " << R << endl;
            break;
         }
      }

      if( polytopeIsValid ) {
         break;
      }

      if( R > 1e16 ) {
         cerr << "Error: could not find valid initial radii" << endl;
         return false;
      }

      R *= initRadiusFactor;
   }
   cout << "Setting initial radii to " << R << endl;

   if( !isEmbeddingConvex() ) {
      cerr << "Warning: initial generalized polyhedral metric is not convex." << endl;
   }

   return true;
}

// Solve the nonlinear equation
//
//    κ(r) - κ* = 0
//
// for r, using Newton's method.  Specifically,
// starting from the current radii r, we iterate
// the expression
//
//    r ← r − J⁻¹|ᵣ(κ(r)−κ*),
//
// where J⁻¹|ᵣ is the Jacobian ∂κ/∂r evaluated at
// the current r.  The values of r are updated
// in-place, and the geometry is updated to reflect
// this change in radii.
//
// Returns true if and only if the solve succeeded.
bool ConvexEmbedder::solveNewton( VertexData<double> kappaStar )
{
   const double tolerance = 1e-4;
   const int maxNewtonIterations = 20;
   size_t n = imesh->nVertices();

   int iteration;
   for( iteration = 0; iteration < maxNewtonIterations; iteration++ )
   {
      // Evaluate u = κ(r)−κ*
      Vector<double> u( n );
      for( Vertex i : imesh->vertices() ) {
         size_t I = i.getIndex();
         u[I] = kappa(i) - kappaStar[I];
      }

      // Check for convergence
      if( u.norm() < tolerance )
      {
         cout << "Newton solver converged to a tolerance of ";
         cout << tolerance << " in ";
         cout << iteration << " iterations" << endl;
         break;
      }

      // Evaluate J⁻¹|ᵣ 
      SparseMatrix<double> J;
      bool failed;
      tie( J, failed ) = buildHessian();

      if( failed )
      {
         cerr << "Hessian has invalid entries; shrinking time step..." << endl;
         return false;
      }

      // Solve J|ᵣ v = u
      Vector<double> v = solveSquare( J, u );

      // Subtract r ← r − τv
      double tau = 1.;
      VertexData<double> r0 = r;
      int maxLineSearchSteps = 32;
      for( int step = 0; step < maxLineSearchSteps; step++ )
      {
         for( Vertex i : imesh->vertices() ) {
            size_t I = i.getIndex();
            r[i] = r0[i] - tau*v[I];
         }
         if( minTriSlack() < 1e-3 ) { // triangles fail to strictly satisfy triangle inequality
            r = r0;
            tau *= .5;
         } else {
            break;
         }
      }
   }

   if( iteration == maxNewtonIterations )
   {
      cout << "Warning: Newton solver failed to ";
      cout << "converge to a tolerance of ";
      cout << tolerance << " after ";
      cout << iteration << " iterations" << endl;
      return false;
   }
   return true; // Newton's method converged
}

// Use Newton's method to find radii with prescribed curvature
void ConvexEmbedder::stepEmbedding()
{
   // Perform gradient descent κₙ₊₁ = κₙ - τ κₙ on the objective
   // ‖κ‖², using a backtracing line search to determine the step
   // size τ, and using Newton's method to find radii rₙ₊₁ corresponding
   // to the next step κₙ₊₁.  If Newton's method fails, we take a
   // smaller step.
   
   VertexData<double> r0 = r; // radii at start of current step
   double tau = options.initialStepSize; // current step size
   const double beta = .5; // backtracking parameter

   for( int step = 0; step < options.maxLineSearchSteps; step++ )
   {
      // Compute next desired curvatures
      VertexData<double> kappaStar( *imesh );
      for( Vertex i : imesh->vertices() )
      {
         kappaStar[i] = (1.-tau)*kappa(i);
      }

      // Attempt to solve for radii r that exhibit
      // the curvatures kappaStar, using Newton's method.
      // (This method updates the radii r in place.)
      bool success = solveNewton( kappaStar );

      // finalize this step if successful
      if( success ) {
         flipToConvex();
         break; // done with line search
      } else { // otherwise, try again
         r = r0; // reset the radii
         tau *= beta; // shrink the time step
      }
   }
}

double ConvexEmbedder::hessianOffDiagonal( Halfedge ij )
{
   Halfedge ji = ij.twin(); 

   return (cot(alpha(ij)) + cot(alpha(ji))) /
          (l(ij) * sin(rho(ij)) * sin(rho(ji)));
}

double ConvexEmbedder::hessianDiagonal( Vertex i )
{
   double sum = 0.;
   for( Halfedge ji : i.incomingHalfedges() ) {
      Vertex j = ji.vertex();
      Halfedge ij = ji.twin();
      if( j == i ) { // self-edge
         sum += l(ji) * (cot(alpha(ji)) + cot(alpha(ij))) / (r[i]*r[i] * sin(rho(ji))*sin(rho(ji))) / 2.;
      } else { // ordinary edge ( i != j )
         sum += cos(phi(ij)) * hessianOffDiagonal( ij );
      }
   }
   return -sum;
}

pair< SparseMatrix<double>, bool >
ConvexEmbedder::buildHessian()
{
   typedef Eigen::Triplet<double> entry;
   vector<entry> entries;
   bool failed = false;

   // Off-diagonal entries.
   // Note that multiple edges may have the same
   // endpoints, which is accounted for by the fact
   // that Eigen will automatically sum contributions
   // from entries with the same row/column indices.
   for( Halfedge ij : imesh->halfedges() )
   {
      size_t I = ij.tipVertex().getIndex();
      size_t J = ij.tailVertex().getIndex();
      if( I != J ) { // self-edges don't correspond to off-diagonals
         double Hij = hessianOffDiagonal( ij );
         entries.push_back( entry( I, J, Hij ));
         if( isinf(Hij) || isnan(Hij) ) {
            failed = true;
            break;
         }
      }
   }

   // Diagonal entries.
   // Unlike off-diagonals, there is only ever one
   // vertex corresponding to each diagonal entry.
   for( Vertex i : imesh->vertices() )
   {
      size_t I = i.getIndex();
      double Hii = hessianDiagonal( i );
      entries.push_back( entry( I, I, Hii ));
      if( isinf(Hii) || isnan(Hii) ) {
         failed = true;
         break;
      }
   }

   size_t n = imesh->nVertices();
   SparseMatrix<double> hessH( n, n );

   if( !failed ) {
      hessH.setFromTriplets( entries.begin(), entries.end() );
   }

   return pair<SparseMatrix<double>,bool>( hessH, failed );
}


double ConvexEmbedder::minTriSlack() {
   double m = numeric_limits<double>::max();

   for( Halfedge ij : imesh->halfedges() )
   {
      Vertex i = ij.vertex();
      Vertex j = ij.next().vertex();

      m = min( m, l(ij) + r[i]  - r[j]  );
      m = min( m, r[i]  + r[j]  - l(ij) );
      m = min( m, r[j]  + l(ij) - r[i]  );
   }

   return m;
}

// Returns the interior angle of the boundary triangle ijk at corner i
double ConvexEmbedder::gamma( Corner ijk )
{
   Halfedge ij = ijk.halfedge();
   Halfedge jk = ij.next();
   Halfedge ki = jk.next();

   return interiorAngle( l(jk), l(ij), l(ki) );
}

double ConvexEmbedder::cosGamma( Corner ijk )
{
   Halfedge ij = ijk.halfedge();
   Halfedge jk = ij.next();
   Halfedge ki = jk.next();

   return interiorAngleCosine( l(jk), l(ij), l(ki) );
}

// Returns the cosine and sine of gamma(ijk)
pair<double,double> ConvexEmbedder::cosSinGamma( Corner ijk )
{
   Halfedge ij = ijk.halfedge();
   Halfedge jk = ij.next();
   Halfedge ki = jk.next();

   double c = max( -1., min( 1., interiorAngleCosine( l(jk), l(ij), l(ki) )));
   double s = sqrt( max( 0., 1. - c*c ));

   return pair<double,double>( c, s );
}


// Returns the angle between boundary edge ij and the radius aj
double ConvexEmbedder::rho( Halfedge ij )
{
   Vertex i = ij.vertex();
   Vertex j = ij.next().vertex();

   return interiorAngle( r[i], l(ij), r[j] );
}

// Returns the cosine of rho(ij)
double ConvexEmbedder::cosRho( Halfedge ij )
{
   Vertex i = ij.vertex();
   Vertex j = ij.next().vertex();

   return interiorAngleCosine( r[i], l(ij), r[j] );
}

// Returns the cosine and sine of rho(ij)
pair<double,double> ConvexEmbedder::cosSinRho( Halfedge ij )
{
   Vertex i = ij.vertex();
   Vertex j = ij.next().vertex();

   double c = max( -1., min( 1., interiorAngleCosine( r[i], l(ij), r[j] )));
   double s = sqrt( max( 0., 1. - c*c ));

   return pair<double,double>( c, s );
}

// Returns the interior angle at vertex a of triangle aij (where a is the apex)
double ConvexEmbedder::phi( Halfedge ij )
{
   Vertex i = ij.vertex();
   Vertex j = ij.next().vertex();

   return interiorAngle( l(ij), r[i], r[j] );
}

// Returns the intrinsic length of edge ij
double ConvexEmbedder::l( Edge ij )
{
   return igeom->edgeLengths[ij];
}

// Returns the intrinsic length of ij.edge()
double ConvexEmbedder::l( Halfedge ij )
{
   return igeom->edgeLengths[ij.edge()];
}

// Returns the dihedral angle around boundary edge ij, on the side containing the halfedge
double ConvexEmbedder::alpha( Halfedge ij )
{
   Halfedge ji = ij.twin();
   Halfedge ki = ij.next().next();
   Corner ijk = ij.corner();

   double cRhoKI = cosRho(ki);

   double cRhoJI, sRhoJI;
   tie( cRhoJI, sRhoJI ) = cosSinRho(ji);

   double cGammaIJK, sGammaIJK;
   tie( cGammaIJK, sGammaIJK ) = cosSinGamma(ijk);

   return acosClamp( (cRhoKI-cRhoJI*cGammaIJK) / (sRhoJI*sGammaIJK) );
}

// Returns the dihedral angle around radius ai, between j and k
double ConvexEmbedder::beta( Corner ijk )
{
   Halfedge ji = ijk.halfedge().twin();
   Halfedge ki = ijk.halfedge().next().next();

   double cGammaIJK = cosGamma( ijk );
   
   double cRhoJI, sRhoJI;
   tie( cRhoJI, sRhoJI ) = cosSinRho(ji);

   double cRhoKI, sRhoKI;
   tie( cRhoKI, sRhoKI ) = cosSinRho(ki);

   return acosClamp( (cGammaIJK-cRhoJI*cRhoKI) / (sRhoJI*sRhoKI) );
}

double ConvexEmbedder::theta( Edge ij )
{
   double sum = 0.;
   for( Halfedge h : ij.adjacentHalfedges() )
   {
      sum += alpha(h);
   }
   return sum;
}

double ConvexEmbedder::omega( Vertex i )
{
   double sum = 0.;
   for( Corner c : i.adjacentCorners() )
   {
      sum += beta(c);
   }
   return sum;
}

double ConvexEmbedder::delta( Vertex i )
{
   double sum = 2.*M_PI;
   for( Corner c : i.adjacentCorners() )
   {
      sum -= gamma(c);
   }
   return sum;
}

double ConvexEmbedder::kappa( Vertex i )
{
   return 2.*M_PI - omega(i);
}

double ConvexEmbedder::cot( double angle )
{
   return -tan( angle + M_PI/2. );
}

double ConvexEmbedder::acosClamp( double angle )
{
   double angleClamp = max( -1., min( 1., angle ));
   return acos( angleClamp );
}

// Returns the cosine of the angle opposite edge a, assuming a,b,c satisfy the triangle inequality.
double ConvexEmbedder::interiorAngleCosine( double a, double b, double c )
{
   return (b*b + c*c - a*a)/(2.*b*c);
}

// Returns the angle opposite edge a, assuming a,b,c satisfy the triangle inequality.
double ConvexEmbedder::interiorAngle( double a, double b, double c )
{
   return acosClamp( interiorAngleCosine( a, b, c ));
}

double ConvexEmbedder::circumradius( Face f )
{
   Halfedge ij = f.halfedge();
   Halfedge jk = ij.next();
   Halfedge ki = jk.next();

   double A = triangleArea( l(ij), l(jk), l(ki) );
   double R = l(ij) * l(jk) * l(ki) / ( 4.*A );

   return R;
}

// Flips any edge with a negative dihedral angle
void ConvexEmbedder::flipToConvex()
{
   const double eps = 1e-3; // tolerance for convexity

   // Enqueue all edges, and keep track
   // of which edges are in the queue
   queue<Edge> Q;
   EdgeData<bool> inQ( *imesh, true );
   for( Edge ij : imesh->edges() ) {
      Q.push( ij );
   }

   // Flip nonconvex edges until none remain
   int nFlips = 0;
   while( !Q.empty() ) {
      Edge ij = Q.front(); Q.pop();
      inQ[ij] = false;

      if( theta(ij) > M_PI + eps ) {
         if( igeom->flipEdgeIfPossible(ij) ) {
            // Enqueue neighbors if not already in queue
            Edge jk = ij.halfedge().next().edge();
            Edge ki = ij.halfedge().next().next().edge();
            Edge il = ij.halfedge().twin().next().edge();
            Edge lj = ij.halfedge().twin().next().next().edge();
            if( jk!=ij  && !inQ[jk] ) { Q.push( jk ); inQ[jk] = true; }
            if( ki!=ij  && !inQ[ki] ) { Q.push( ki ); inQ[ki] = true; }
            if( il!=ij  && !inQ[il] ) { Q.push( il ); inQ[il] = true; }
            if( lj!=ij  && !inQ[lj] ) { Q.push( lj ); inQ[lj] = true; }
            nFlips++;
         } else {
            cerr << "Warning: failed to flip nonconvex edge!" << endl;
         }
      }
   }

   cout << "Flipped " << nFlips << " edges." << endl;
}

void ConvexEmbedder::refreshVertexCoordinates() {

   // Compute coordinates at each triangle corner by incrementally
   // laying out tetrahedra defined by the radii.  Note that if the
   // generalized polyhedral metric has nonzero curvature along radial
   // edges, these tetrahedra will not all agree in three-dimensional space.
   localLayout = CornerData<Vector3>( *imesh );
   CornerData<Vector3>& x( localLayout ); // shorthand for "localLayout"

   // Embed the tetrahedron corresponding to the first triangle
   Face ijk = imesh->face(0);
   Halfedge ij = ijk.halfedge();
   Halfedge jk = ij.next();
   Halfedge ki = jk.next();
   Corner i = ij.corner();
   Corner j = jk.corner();
   Corner k = ki.corner();
   double Lai = r[i.vertex()];
   double Laj = r[j.vertex()];
   double Lak = r[k.vertex()];
   double cosThetaI = interiorAngleCosine( l(jk), l(ki), l(ij) );
   double sinThetaI = sqrt( 1. - cosThetaI*cosThetaI );
   x[i] = Vector3{ 0., 0., 0. };
   x[j] = Vector3{ l(ij), 0., 0. };
   x[k] = Vector3{ l(ki)*cosThetaI, l(ki)*sinThetaI, 0. };
   apex = tetFourthPoint( x[i], x[j], x[k], Lai, Laj, Lak );

   // Perform breadth-first traversal to incrementally embed neighboring tets
   queue<Halfedge> Q;
   FaceData<bool> visited( *imesh, false ); // keep track of which tets have been visited
   visited[ijk] = true;
   if( !visited[ij.twin().face() ] ) { Q.push( ij.twin() ); visited[ij.twin().face()] = true; }
   if( !visited[jk.twin().face() ] ) { Q.push( jk.twin() ); visited[jk.twin().face()] = true; }
   if( !visited[ki.twin().face() ] ) { Q.push( ki.twin() ); visited[ki.twin().face()] = true; }
   while( !Q.empty() )
   {
      ij = Q.front(); Q.pop();
      jk = ij.next();
      ki = jk.next();
      i = ij.corner();
      j = jk.corner();
      k = ki.corner();
      ijk = ij.face();

      x[i] = x[ij.twin().next().corner()];
      x[j] = x[ij.twin().corner()];
      x[k] = tetFourthPoint( x[i], apex, x[j], l(ki), r[k.vertex()], l(jk) );

      if( !visited[jk.twin().face()] ) { Q.push( jk.twin() ); visited[jk.twin().face()] = true; }
      if( !visited[ki.twin().face()] ) { Q.push( ki.twin() ); visited[ki.twin().face()] = true; }
   }

   // Average corner values to final vertex values
   embedding = VertexData<Vector3>( *imesh );
   for( Vertex i : imesh->vertices() )
   {
      embedding[i] = Vector3{0.,0.,0.};
      for( Corner ijk : i.adjacentCorners() )
      {
         embedding[i] += localLayout[ijk];
      }
      embedding[i] /= i.degree();
   }
}

} // namespace surface
} // namespace geometrycentral
