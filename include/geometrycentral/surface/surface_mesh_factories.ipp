#pragma once


namespace geometrycentral {
namespace surface {

template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<ManifoldSurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeManifoldSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat) {

  std::unique_ptr<ManifoldSurfaceMesh> mesh(new ManifoldSurfaceMesh(fMat));
  std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh, vMat));

  return std::make_tuple(std::move(mesh), std::move(geometry));
}

template <typename Scalar_V, typename Scalar_F>
std::tuple<std::unique_ptr<SurfaceMesh>, std::unique_ptr<VertexPositionGeometry>>
makeSurfaceMeshAndGeometry(const Eigen::MatrixBase<Scalar_V>& vMat, const Eigen::MatrixBase<Scalar_F>& fMat) {

  std::unique_ptr<SurfaceMesh> mesh(new SurfaceMesh(fMat));
  std::unique_ptr<VertexPositionGeometry> geometry(new VertexPositionGeometry(*mesh, vMat));

  return std::make_tuple(std::move(mesh), std::move(geometry));
}

} // namespace surface
} // namespace geometrycentral
