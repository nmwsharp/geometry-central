`geometrycentral::Vector3` is the basic 3D vector type in geometry central. There are many like it, but this one is ours.

`#!cpp #include "geometrycentral/utilities/vector3.h"`

### Construction 

`Vector3` is a POD type, so you should use brace-initialization sytax:

```cpp
#include "geometrycentral/utilities/vector3.h
using namespace geometrycentral;

Vector3 myVec{3.8, 2.9, 1.1}; //create
myVec = Vector3{1.1, 2.2, 3.3}; // reassign
```

Factory methods can construct a few common values:

??? func "`#!cpp static Vector3 Vector3::zero()`"

    Returns the zero vector

??? func "`#!cpp static Vector3 Vector3::constant(double c)`"

    Returns a vector with all components set to $c$


??? func "`#!cpp static Vector3 Vector3::infinity()`"

    Returns the infinite vector $(\infty, \infty, \infty)$.

??? func "`#!cpp static Vector3 Vector3::undefined()`"

    Returns the undefined vector `(NaN, NaN, NaN)`.


### Access

The three elements of the vector can be accessed as `vec.x` and `vec.y` and `vec.z`.

Alternately, the elements can be indexed as `vec[0]` and `vec[1]` and `vec[2]`.


### Conversion

??? func "`#!cpp Vector3::operator<<()`"

    `Vector3` can be serialized.

    ```cpp
    Vector3 v{1.2, 3.4, 5.6};
    std::cout << v << std::endl;
    // prints something like: <1.2, 3.4, 5.6>
    ```

### Arithmetic

Vector3 supports the element-wise addition, subtraction, and scalar multiplication you would probably expect.


### Member operations

These methods _do not_ change the underlying `Vector3`, but return a new `Vector3`.
```cpp
Vector3 vec{1., 2., 3.};
vec.normalize();        // does nothing
vec = vec.normalize();  // much better
```

??? func "`#!cpp Vector3 Vector3::normalize()`"

    Returns a unit-norm vector pointing in the same direction. If the input is the zero vector, the result will contain NaNs.

??? func "`#!cpp Vector3 Vector3::normalizeCutoff(double mag = 0.)`"

    Returns a unit-norm vector with the same direction. If the input has magnitude less `<= mag`, the vector is unchanged.

??? func "`#!cpp Vector3 Vector3::unit()`"

    Alias for `normalize()`.

??? func "`#!cpp Vector3 Vector3::rotateAround(Vector3 axis, double theta)`"

    Rotate the vector by angle $\theta$ around `axis` in the right-handed direction. `axis` need not be a unit vector.

??? func "`#!cpp Vector3 Vector3::removeComponent(Vector3 unitDir)`"

    Removes any component of this vector in the direction `unitDir`, making the result orthogonal to `unitDir`. As the name suggests, `unitDir` must be a unit vector.

??? func "`#!cpp double Vector3::norm()`"

    Returns the magnitude of the vector.

    Also available as `norm(v)`.


??? func "`#!cpp double Vector3::norm2()`"

    Returns the squared magnitude of the vector.

    Also available as `norm2(v)`.


### Function operations

These operations do not change the vector on which they are called.

??? func "`#!cpp double norm(Vector3 v)`"

    Returns the magnitude of the vector.

    Also available as `v.norm()`.


??? func "`#!cpp double norm2(Vector3 v)`"

    Returns the squared magnitude of the vector.

    Also available as `v.norm2()`.

??? func "`#!cpp Vector3 normalize(Vector3 v)`"

    Returns normalized copy of the vector.

??? func "`#!cpp Vector3 normalizeCutoff(Vector3 v, double mag = 0.)`"

    Returns a normalized copy of the vector. If the input has magnitude less `<= mag`, the vector is unchanged.

??? func "`#!cpp Vector3 unit(Vector3 v)`"

    Alias for `normalize(v)`.

??? func "`#!cpp double dot(Vector3 u, Vector3 v)`"

    Returns the dot product between two vectors.


??? func "`#!cpp double sum(Vector3 u)`"

    Returns the sum of the coordinates of a vector


??? func "`#!cpp Vector3 cross(Vector3 u, Vector3 v)`"

    Returns the cross product between two vectors.


??? func "`#!cpp double angle(Vector3 u, Vector3 v)`"

    Returns the angle between two not-necessarily-unit vectors. Output in the range $[0, \pi]$.


??? func "`#!cpp double angleInPlane(Vector3 u, Vector3 v, Vector3 normal)`"

    Returns the signed angle between two not-necessarily-unit vectors, measured in the plane defined by `normal` (which need not be a unit vector). Output is in the range $[-\pi, \pi]$, as in `atan2`.


??? func "`#!cpp Vector3 clamp(Vector3 val, Vector3 low, Vector3 high)`"

    Returns returns a a vector where each component has been clamped to be between the corresponding compnents of `low` and `high`.

??? func "`#!cpp Vector3 componentwiseMin(Vector3 u, Vector3 v)`"

    Returns a new vector, each component of which is the minimum of that component in `u` and `v`.


??? func "`#!cpp Vector3 componentwiseMax(Vector3 u, Vector3 v)`"

    Returns a new vector, each component of which is the maximum of that component in `u` and `v`.


### Properties

??? func "`#!cpp bool isfinite(Vector3 u)`"

    Returns true if all of the components of the vector are finite.

    Note: this function is intentionally not camel-cased out of solidarity with `std::isfinite()`.

    Also available as `u.isFinite()`.


??? func "`#!cpp bool isDefined(Vector3 u)`"

    Returns true if all of the components of the vector are not NaN.
    
    Also available as `u.isDefined()`.
