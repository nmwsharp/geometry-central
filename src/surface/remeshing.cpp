#include "geometrycentral/surface/remeshing.h"

namespace geometrycentral {
namespace surface {

using std::queue;

Vector3 vertexNormal(VertexPositionGeometry& geometry, Vertex v)
{
    Vector3 norm = Vector3::zero();
    for(Corner c : v.adjacentCorners()){
        norm += geometry.cornerAngle(c) * geometry.faceNormal(c.face());
    }
    return normalize(norm);
}

inline Vector3 projectToPlane(Vector3 v, Vector3 norm)
{
    return v - norm * dot(norm, v);
}

inline Vector3 edgeMidpoint(SurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e)
{
    Vector3 endPos1 = geometry.inputVertexPositions[e.halfedge().tailVertex()];
    Vector3 endPos2 = geometry.inputVertexPositions[e.halfedge().tipVertex()];
    return (endPos1+endPos2)/2;
}

Vector3 findCircumcenter(Vector3 p1, Vector3 p2, Vector3 p3)
{
    // barycentric coordinates of circumcenter
    double a = (p3 - p2).norm();
    double b = (p3 - p1).norm();
    double c = (p2 - p1).norm();
    double a2 = a * a;
    double b2 = b * b;
    double c2 = c * c;
    Vector3 O{a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2), c2 * (a2 + b2 - c2)};
    // normalize to sum of 1
    O /= O[0] + O[1] + O[2];
    // change back to space
    return O[0] * p1 + O[1] * p2 + O[2] * p3;
}

Vector3 findCircumcenter(VertexPositionGeometry& geometry, Face f)
{
    // retrieve the face's vertices
    int index = 0;
    Vector3 p[3];
    for (Vertex v0 : f.adjacentVertices())
    {
        p[index] = geometry.inputVertexPositions[v0];
        index++;
    }
    return findCircumcenter(p[0], p[1], p[2]);
}

bool isDelaunay(VertexPositionGeometry& geometry, Edge e)
{
    float angle1 = geometry.cornerAngle(e.halfedge().next().next().corner());
    float angle2 = geometry.cornerAngle(e.halfedge().twin().next().next().corner());
    return angle1 + angle2 <= PI;
}

inline double diamondAngle(Vector3 a, Vector3 b, Vector3 c, Vector3 d) // dihedral angle at edge a-b
{
    Vector3 n1 = cross(b-a, c-a);
    Vector3 n2 = cross(b-d, a-d);
    return PI-angle(n1, n2);
}

inline bool checkFoldover(Vector3 a, Vector3 b, Vector3 c, Vector3 x, double angle)
{
    return diamondAngle(a, b, c, x) < angle;
}

bool shouldCollapse(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e)
{
    std::vector<Halfedge> toCheck;
    Vertex v1 = e.halfedge().vertex();
    Vertex v2 = e.halfedge().twin().vertex();
    Vector3 midpoint = edgeMidpoint(mesh, geometry, e);
    // find (halfedge) link around the edge, starting with those surrounding v1
    Halfedge he = v1.halfedge();
    Halfedge st = he;
    do{
        he = he.next();
        if(he.vertex() != v2 && he.next().vertex() != v2){
            toCheck.push_back(he);
        }
        he = he.next().twin();
    }
    while(he != st);
    // link around v2
    he = v2.halfedge();
    st = he;
    do{
        he = he.next();
        if(he.vertex() != v1 && he.next().vertex() != v1){
            toCheck.push_back(he);
        }
        he = he.next().twin();
    }
    while(he != st);
    
    // see if the point that would form after a collapse would cause a major foldover with surrounding edges
    for(Halfedge he0 : toCheck){
        Halfedge heT = he0.twin();
        Vertex v1 = heT.vertex();
        Vertex v2 = heT.next().vertex();
        Vertex v3 = heT.next().next().vertex();
        Vector3 a = geometry.inputVertexPositions[v1];
        Vector3 b = geometry.inputVertexPositions[v2];
        Vector3 c = geometry.inputVertexPositions[v3];
        if(checkFoldover(a, b, c, midpoint, 2)){
            return false;
        }
    }
    return true;
}

double getSmoothMeanCurvature(VertexPositionGeometry& geometry, Vertex v)
{
    double A = geometry.vertexDualArea(v);
    double S = geometry.vertexMeanCurvature(v);
    double K = S / A;
    return K;
}

// flatLength: specifies how long the target edge length should be in flat regions
// epsilon: controls how much variation in target length occurs due to curvature

double findMeanTargetL(SurfaceMesh& mesh, VertexPositionGeometry& geometry, Edge e, double flatLength, double epsilon)
{
    double averageH = 0;
    for (Vertex v : e.adjacentVertices()) {
        averageH += getSmoothMeanCurvature(geometry, v);
    }
    averageH /= 2;
    double L = flatLength * epsilon / (fabs(averageH) + epsilon);
    return L;
}



void fixDelaunay(SurfaceMesh& mesh, VertexPositionGeometry& geometry)
{
    // queue of edges to check if Delaunay
    queue<Edge> toCheck;
    // true if edge is currently in toCheck
    EdgeData<bool> inQueue(mesh);
    // start with all edges
    for (Edge e : mesh.edges())
    {
        toCheck.push(e);
        inQueue[e] = true;
    }
    // counter and limit for number of flips
    int flipMax = 100 * mesh.nVertices();
    int flipCnt = 0;
    while (!toCheck.empty() && flipCnt < flipMax)
    {
        Edge e = toCheck.front();
        toCheck.pop();
        inQueue[e] = false;
        // if not Delaunay, flip edge and enqueue the surrounding "diamond" edges (if not already)
        if (!e.isBoundary() && !isDelaunay(geometry, e))
        {
            flipCnt++;
            Halfedge he = e.halfedge();
            Halfedge he1 = he.next();
            Halfedge he2 = he1.next();
            Halfedge he3 = he.twin().next();
            Halfedge he4 = he3.next();

            if (!inQueue[he1.edge()])
            {
                toCheck.push(he1.edge());
                inQueue[he1.edge()] = true;
            }
            if (!inQueue[he2.edge()])
            {
                toCheck.push(he2.edge());
                inQueue[he2.edge()] = true;
            }
            if (!inQueue[he3.edge()])
            {
                toCheck.push(he3.edge());
                inQueue[he3.edge()] = true;
            }
            if (!inQueue[he4.edge()])
            {
                toCheck.push(he4.edge());
                inQueue[he4.edge()] = true;
            }
            mesh.flip(e);
        }
    }
}

void smoothByLaplacian(SurfaceMesh& mesh, VertexPositionGeometry& geometry)
{
    // smoothed vertex positions
    VertexData<Vector3> newVertexPosition(mesh);
    for (Vertex v : mesh.vertices())
    {
        if(v.isBoundary())
        {
            newVertexPosition[v] = geometry.inputVertexPositions[v];
        }
        else
        {
            // calculate average of surrounding vertices
            newVertexPosition[v] = Vector3::zero();
            for (Vertex j : v.adjacentVertices())
            {
                newVertexPosition[v] += geometry.inputVertexPositions[j];
            }
            newVertexPosition[v] /= v.degree();
            // and project the average to the tangent plane
            Vector3 updateDirection = newVertexPosition[v] - geometry.inputVertexPositions[v];
            updateDirection = projectToPlane(updateDirection, vertexNormal(geometry, v));
            
            newVertexPosition[v] = geometry.inputVertexPositions[v] + 1*updateDirection;
        }
    }
    // update final vertices
    for (Vertex v : mesh.vertices())
    {
        geometry.inputVertexPositions[v] = newVertexPosition[v];
    }
}

void smoothByCircumcenter(SurfaceMesh& mesh, VertexPositionGeometry& geometry)
{
    geometry.requireFaceAreas();
    // smoothed vertex positions
    VertexData<Vector3> newVertexPosition(mesh);
    for (Vertex v : mesh.vertices())
    {
        newVertexPosition[v] = geometry.inputVertexPositions[v]; // default
        if(!v.isBoundary())
        {
            Vector3 updateDirection = Vector3::zero();
            for (Face f : v.adjacentFaces())
            {
                // add the circumcenter weighted by face area to the update direction
                Vector3 circum = findCircumcenter(geometry, f);
                updateDirection += geometry.faceArea(f) * (circum - geometry.inputVertexPositions[v]);
            }
            updateDirection /= (3 * geometry.vertexDualArea(v));
            // project update direction to tangent plane
            updateDirection = projectToPlane(updateDirection, vertexNormal(geometry, v));
            Vector3 newPos = geometry.inputVertexPositions[v] + .5 * updateDirection;
            newVertexPosition[v] = newPos;
        }
    }
    // update final vertices
    for (Vertex v : mesh.vertices())
    {
        geometry.inputVertexPositions[v] = newVertexPosition[v];
    }
}


bool adjustEdgeLengths(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geometry, double flatLength, double epsilon, double minLength, bool curvatureAdaptive)
{
    bool didSplitOrCollapse = false;
    // queues of edges to CHECK to change
    std::vector<Edge> toSplit;
    std::vector<Edge> toCollapse;
    
    for(Edge e : mesh.edges())
    {
        toSplit.push_back(e);
    }
    
    // actually splitting
    while(!toSplit.empty())
    {
        Edge e = toSplit.back();
        toSplit.pop_back();
        double length_e = geometry.edgeLength(e);
        double threshold = (curvatureAdaptive) ? findMeanTargetL(mesh, geometry, e, flatLength, epsilon) : flatLength;
        if(length_e > minLength && length_e > threshold * 1.5)
        {
            Vector3 newPos = edgeMidpoint(mesh, geometry, e);
            Halfedge he = mesh.splitEdgeTriangular(e);
            didSplitOrCollapse = true;
            Vertex newV = he.vertex();
            geometry.inputVertexPositions[newV] = newPos;
        }
        else
        {
            toCollapse.push_back(e);
        }                
        
    }
    // actually collapsing
    while(!toCollapse.empty())
    {
        Edge e = toCollapse.back();
        toCollapse.pop_back();
        if(e.halfedge().next().getIndex() != INVALID_IND) // make sure it exists
        {
            double threshold = (curvatureAdaptive) ? findMeanTargetL(mesh, geometry, e, flatLength, epsilon) : flatLength;
            if(geometry.edgeLength(e) < threshold * 0.5)
            {
                Vector3 newPos = edgeMidpoint(mesh, geometry, e);
                if(shouldCollapse(mesh, geometry, e)) {
                    Vertex v = mesh.collapseEdgeTriangular(e);
                    didSplitOrCollapse = true;
                    if (v != Vertex()) {
                        if(!v.isBoundary()) {
                            geometry.inputVertexPositions[v] = newPos;
                        }
                    }
                }
            }
        }
    }
    return didSplitOrCollapse;
}

} // namespace surface
} // namespace geometrycentral
