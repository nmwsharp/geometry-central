
## Eigen

Generally, geometry central uses [Eigen](http://eigen.tuxfamily.org) for all matrix types. Though we build additional solvers and utilities on top of Eigen.  See [the Eigen section of dependencies](../../build/dependencies/#eigen) for instructions about getting Eigen and integrating with existing build systems.

Note that the [Vector2](../../utilities/vector2/) and [Vector3](../../utilities/vector3) low-dimensional scalar types are entirely separate from these high-dimensional linear algebra types; `Vector2` and `Vector3` do not use Eigen, and they cannot participate in arithmetic expressions with Eigen types.

Two typedefs are used extensively throughout geometry central to make the default Eigen types slightly less verbose. Both are defined in `linear_algebra_utilities.h`.

`#!cpp #include "geometrycentral/numerical/linear_algebra_utilities.h"`

??? func "`#!cpp Vector<T>`"

    A templated vector typedef, to Eigen's vector type.
    ```cpp
    template <typename T>
    using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    ```

    Use like `Vector<double>` or `Vector<bool>`.
    

??? func "`#!cpp SparseMatrix<T>`"

    A templated sparse matrix typedef, to Eigen's sparse matrix type.
    ```cpp
    template <typename T>
    using SparseMatrix = Eigen::SparseMatrix<T>;
    ```

    Use like `SparseMatrix<double>` or `SparseMatrix<int>`.



#### Gotchas

Be wary, Eigen's alignment rules make it efficient, but also impose requirements which can lead to hard-to-debug memory errors. A few particularly common pitfalls are:

- Avoid intermingling different versions of Eigen in the same program. Suppose some part of you codebase uses one version of Eigen, and a dependency uses a different version. Linking a function which returns an Eigen vector between these versions can lead to segfaults, because different alignment policies were used. 

- Similar to the previous, linking Eigen programs compiled with different preprocessor directives and optimization flags can yield binary incompatibility. Be sure that all parts of your codebase using Eigen receive the same build options.

- Fixed-sized Eigen types (like `Eigen::Vector4d`, but **not** our `Vector<T>` or `SparseMatrix<T>`) *may not* be passed by value to functions. The same applies transitively to classes which have Eigen types as members. See this (opinionated) [note](https://eigen.tuxfamily.org/dox-devel/group__TopicPassingByValue.html).
