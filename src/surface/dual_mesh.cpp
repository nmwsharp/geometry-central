#include "geometrycentral/surface/dual_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;


std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
dual_mesh_nokeep(ManifoldSurfaceMesh& mesh,
                 VertexPositionGeometry& geometry)
{
    std::vector<Vector3> Barycs;
    Barycs.resize(mesh.nFaces());
    for (Face f : mesh.faces())
    {
        Barycs[f.getIndex()] = { 0.0, 0.0, 0.0 };
        for (Vertex vf : f.adjacentVertices())
            Barycs[f.getIndex()] += geometry.vertexPositions[vf];
        Barycs[f.getIndex()] /= f.degree();
    }

    std::vector<std::vector<size_t>> VFaces;
    VFaces.reserve(mesh.nVertices());
    for (Vertex v : mesh.vertices())
    {
        std::vector<size_t> Poly;
        Poly.reserve(v.faceDegree() + (v.isBoundary() ? 3 : 0));
        for (Face fv : v.adjacentFaces())
            Poly.emplace_back(fv.getIndex());
        if (v.isBoundary())
            continue;
        std::reverse(Poly.begin(), Poly.end());
        VFaces.emplace_back(Poly);
    }

    return makeManifoldSurfaceMeshAndGeometry(VFaces, Barycs);
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
dual_mesh_keep(ManifoldSurfaceMesh& mesh,
               VertexPositionGeometry& geometry)
{
    size_t ExtraVerts = 0;
    if (mesh.nBoundaryLoops() > 0)
    {
        for (BoundaryLoop bl : mesh.boundaryLoops())
        {
            for (Vertex v : bl.adjacentVertices())
                ExtraVerts++;
            for (Edge e : bl.adjacentEdges())
                ExtraVerts++;
        }
    }

    std::vector<Vector3> Barycs;
    Barycs.reserve(mesh.nFaces() + ExtraVerts);
    Barycs.resize(mesh.nFaces());
    for (Face f : mesh.faces())
    {
        Barycs[f.getIndex()] = { 0.0, 0.0, 0.0 };
        for (Vertex vf : f.adjacentVertices())
            Barycs[f.getIndex()] += geometry.vertexPositions[vf];
        Barycs[f.getIndex()] /= f.degree();
    }
    size_t CurVLen = Barycs.size();

    std::vector<int> XVertsMap;
    XVertsMap.resize(mesh.nVertices(), -1);
    std::vector<std::vector<size_t>> VFaces;
    VFaces.resize(mesh.nVertices());
    for (Vertex v : mesh.vertices())
    {
        auto& Poly = VFaces[v.getIndex()];
        Poly.reserve(v.faceDegree() + (v.isBoundary() ? 3 : 0));
        for (Face fv : v.adjacentFaces())
            Poly.emplace_back(fv.getIndex());
        if (v.isBoundary())
        {
            int Offset = 0;
            for (Face fv : v.adjacentFaces())
            {
                if (fv == v.halfedge().face())
                  break;
                Offset++;
            }
            std::rotate(Poly.begin(), Poly.begin() + Offset + 1, Poly.end());
            if (XVertsMap[v.getIndex()] == -1)
            {
                Barycs.emplace_back(0.5 * (geometry.vertexPositions[v] + geometry.vertexPositions[v.halfedge().next().vertex()]));
                Poly.push_back(CurVLen);
                XVertsMap[v.getIndex()] = CurVLen;
                CurVLen++;
            }
            else
                Poly.push_back(XVertsMap[v.getIndex()]);
            Barycs.emplace_back(geometry.vertexPositions[v]);
            Poly.push_back(CurVLen);
            CurVLen++;
            int fvid = Poly[0];
            Face fv = mesh.face(fvid);
            for (Halfedge he : fv.adjacentHalfedges())
            {
                if (he.next().vertex() != v)
                    continue;
                Vertex ov = he.vertex();
                if (XVertsMap[ov.getIndex()] == -1)
                {
                    Barycs.emplace_back(0.5 * (geometry.vertexPositions[v] + geometry.vertexPositions[ov]));
                    Poly.push_back(CurVLen);
                    XVertsMap[ov.getIndex()] = CurVLen;
                    CurVLen++;
                }
                else
                    Poly.push_back(XVertsMap[ov.getIndex()]);
            }
        }
        std::reverse(Poly.begin(), Poly.end());
    }

    return makeManifoldSurfaceMeshAndGeometry(VFaces, Barycs);
}


























std::unique_ptr<ManifoldSurfaceMesh>
dual_mesh_nokeep(ManifoldSurfaceMesh& mesh)
{
    std::vector<std::vector<size_t>> VFaces;
    VFaces.reserve(mesh.nVertices());
    for (Vertex v : mesh.vertices())
    {
        std::vector<size_t> Poly;
        Poly.reserve(v.faceDegree() + (v.isBoundary() ? 3 : 0));
        for (Face fv : v.adjacentFaces())
            Poly.emplace_back(fv.getIndex());
        if (v.isBoundary())
            continue;
        std::reverse(Poly.begin(), Poly.end());
        VFaces.emplace_back(Poly);
    }

    std::unique_ptr<ManifoldSurfaceMesh> res(new ManifoldSurfaceMesh(VFaces));
    return std::move(res);
}

std::unique_ptr<ManifoldSurfaceMesh>
dual_mesh_keep(ManifoldSurfaceMesh& mesh)
{
    size_t ExtraVerts = 0;
    if (mesh.nBoundaryLoops() > 0)
    {
        for (BoundaryLoop bl : mesh.boundaryLoops())
        {
            for (Vertex v : bl.adjacentVertices())
                ExtraVerts++;
            for (Edge e : bl.adjacentEdges())
                ExtraVerts++;
        }
    }
    
    size_t CurVLen = mesh.nFaces();

    std::vector<int> XVertsMap;
    XVertsMap.resize(mesh.nVertices(), -1);
    std::vector<std::vector<size_t>> VFaces;
    VFaces.resize(mesh.nVertices());
    for (Vertex v : mesh.vertices())
    {
        auto& Poly = VFaces[v.getIndex()];
        Poly.reserve(v.faceDegree() + (v.isBoundary() ? 3 : 0));
        for (Face fv : v.adjacentFaces())
            Poly.emplace_back(fv.getIndex());
        if (v.isBoundary())
        {
            int Offset = 0;
            for (Face fv : v.adjacentFaces())
            {
                if (fv == v.halfedge().face())
                  break;
                Offset++;
            }
            std::rotate(Poly.begin(), Poly.begin() + Offset + 1, Poly.end());
            if (XVertsMap[v.getIndex()] == -1)
            {
                Poly.push_back(CurVLen);
                XVertsMap[v.getIndex()] = CurVLen;
                CurVLen++;
            }
            else
                Poly.push_back(XVertsMap[v.getIndex()]);
            Poly.push_back(CurVLen);
            CurVLen++;
            int fvid = Poly[0];
            Face fv = mesh.face(fvid);
            for (Halfedge he : fv.adjacentHalfedges())
            {
                if (he.next().vertex() != v)
                    continue;
                Vertex ov = he.vertex();
                if (XVertsMap[ov.getIndex()] == -1)
                {
                    Poly.push_back(CurVLen);
                    XVertsMap[ov.getIndex()] = CurVLen;
                    CurVLen++;
                }
                else
                    Poly.push_back(XVertsMap[ov.getIndex()]);
            }
        }
        std::reverse(Poly.begin(), Poly.end());
    }

    std::unique_ptr<ManifoldSurfaceMesh> res(new ManifoldSurfaceMesh(VFaces));
    return std::move(res);
}





















std::unique_ptr<ManifoldSurfaceMesh>
geometrycentral::surface::dual_mesh(ManifoldSurfaceMesh& mesh,
                                   bool keepBoundaries)
{
    if (keepBoundaries)
        return dual_mesh_keep(mesh);
    else
        return dual_mesh_nokeep(mesh);
}

std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
geometrycentral::surface::dual_mesh(ManifoldSurfaceMesh& mesh,
                                   VertexPositionGeometry& geometry,
                                   bool keepBoundaries)
{
    if (keepBoundaries)
        return dual_mesh_keep(mesh, geometry);
    else
        return dual_mesh_nokeep(mesh, geometry);
}