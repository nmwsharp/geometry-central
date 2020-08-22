Helper functions to interface with Eigen data in interesting ways.

`#!cpp #include "geometrycentral/utilities/eigen_interop_helpers.h"`

## Expanding a homogenous POD type into an `Eigen::Matrix`

??? func "`#!cpp Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, k, Options>, Alignment> EigenMap(MeshData<E,O> &data)`"

    Given `MeshData` storing a POD type `O` like `Vector3` which can be meaningfully decomposed into a set of `k` other types `T`, this function maps the memory of the vector of `O` as a `Matrix<T, Dynamic, k>`. E.g., `MeshData<Vertex, Vector3>` (N x 1) -> `Map<Matrix<double, Dynamic, 3>>` (N x 3). The template `Options` allows you to specify is the underlying data should be bound in RowMajor or ColMajor order.

    ```cpp
    VertexData<Vector3> pos;
    // auto resolves to Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>, Alignment>
    auto mapped_position_matrix = EigenMap<double, 3, Eigen::RowMajor>(pos);
    ```

    Note that it is unnecessary to specify `Eigen::RowMajor` in this example since it is the default template option. Furthermore template options `E` and `O` corresponding to the Mesh Element and Data Type need not be explicitly specified since the compiler can infer it for us.

    Most important, modifications to `mapped_position_matrix` will be reflected in `pos` since they share the same memory!

??? func "`#!cpp Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, k, Options>, Alignment> EigenMap(const MeshData<E,O> &data)`"

    Const version of the above.

## Flattening a homogeneous POD type into an `Eigen::Vector`

??? func "`#!cpp Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1, Options>, Alignment> FlattenedEigenMap(MeshData<E,O> &data)`"

    Given `MeshData` storing a POD type `O` like `Vector3` which can be meaningfully decomposed into a set of `k` other types `T`, this function maps the memory of the vector of `O` as a vector of `T` (`Matrix<T, Dynamic, 1>`). E.g., `MeshData<Vertex, Vector3>` (N x 1) -> `Map<Matrix<double, Dynamic, 1>>` (3N x 1). The template `Options` allows you to specify is the underlying data should be bound in RowMajor or ColMajor order.

    ```cpp
    VertexData<Vector3> pos;
    // auto resolves to Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>, Alignment>
    auto mapped_flat_position_matrix = FlattenedEigenMap<double, 3>(pos);
    ```

    Note that the return defaults to `ColMajor` order. The template options `E` and `O` corresponding to the Mesh Element and Data Type need not be explicitly specified since the compiler can infer it for us.

    Most important, modifications to `mapped_flat_position_matrix` will be reflected in `pos` since they share the same memory! 

??? func "`#!cpp Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1, Options>, Alignment> FlattenedEigenMap(const MeshData<E,O> &data)`"

    Const version of the above.
