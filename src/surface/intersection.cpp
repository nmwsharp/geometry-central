#include "geometrycentral/surface/intersection.h"
#include "geometrycentral/utilities/elementary_geometry.h"

#include <queue>
#include <tuple>


namespace geometrycentral {
namespace surface {

SurfaceIntersectionResult selfIntersections( VertexPositionGeometry& geometry ) {

   SurfaceMesh& mesh = geometry.mesh;
   SurfaceIntersectionResult intersections;
   intersections.hasIntersections = false;
   size_t n = 0; // number of intersection points

   // vertices and vertex locations for the current triangle pair
   Vertex u[3];
   Vertex v[3];
   Vector3 p[3];
   Vector3 q[3];

   // iterate over all face pairs f,g
   for( Face f : mesh.faces() ) {

      // get vertices ui of f
      int i = 0;
      for( Vertex ui : f.adjacentVertices() ) {
         u[i] = ui;
         i++;
      }

      for( Face g : mesh.faces() ) {

         // get vertices uj of g
         int j = 0;
         for( Vertex vj : g.adjacentVertices() ) {
            v[j] = vj;
            j++;
         }

         // skip triangles that share vertices
         if( u[0] == v[0] || u[0] == v[1] || u[0] == v[2] ||
             u[1] == v[0] || u[1] == v[1] || u[1] == v[2] ||
             u[2] == v[0] || u[2] == v[1] || u[2] == v[2] ) {
            break;
         }

         // get vertex locations
         for( int i = 0; i < 3; i++ ) {
            p[i] = geometry.vertexPositions[u[i]];
            q[i] = geometry.vertexPositions[v[i]];
         }

         // check for and compute intersection
         TriTriIntersectionResult3D r;
         r = triTriIntersection( p[0], p[1], p[2],
               q[0], q[1], q[2] );

         // add to list of all intersections
         if( r.intersect ) {
            intersections.hasIntersections = true;
            intersections.points.push_back( r.xA );
            intersections.points.push_back( r.xB );
            intersections.edges.push_back( { n, n+1 } );
            n += 2;
         }
      }
   }

   return intersections;
}

} // namespace surface
} // namespace geometrycentral
