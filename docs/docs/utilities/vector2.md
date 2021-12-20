`geometrycentral::Vector2` is the basic 2D vector type in geometry central. Like a good turkey sandwich, it aims to be unsurprising yet satisfying.

Of particular interest, `Vector2` is also used to encode 2D rotations, by supporting multiplication as a complex number. See the [rotations section](#rotations-and-complex-multiplication).

`#!cpp #include "geometrycentral/utilities/vector2.h"`

### Construction 

`Vector2` is a POD type, so you should use brace-initialization sytax:

```cpp
#include "geometrycentral/utilities/vector2.h
using namespace geometrycentral;

Vector2 myVec{3.8, 2.9}; //create
myVec = Vector2{1.1, 2.2}; // reassign
```

Factory methods can construct a few common values:

??? func "`#!cpp static Vector2 Vector2::zero()`"

    Returns the zero vector

??? func "`#!cpp static Vector2 Vector2::constant(double c)`"

    Returns a vector with all components set to $c$

??? func "`#!cpp static Vector2 Vector2::infinity()`"

    Returns the infinite vector $(\infty, \infty)$.

??? func "`#!cpp static Vector2 Vector2::undefined()`"

    Returns the undefined vector `(NaN, NaN)`.

And serve as constructors:

??? func "`#!cpp static Vector2 Vector2::fromAngle(double theta)`"

    Returns the vector $(\cos(\theta), \sin(\theta))$.

??? func "`#!cpp static Vector2 Vector2::fromComplex(std::complex<double> c)`"

    Converts a `std::complex<double>` to a `Vector2`.

### Access

The two elements of the vector can be accessed as `vec.x` and `vec.y`.

Alternately, the two elements can be indexed as `vec[0]` and `vec[1]`.


### Conversion

??? func "`#!cpp Vector2::operator std::complex<double>()`"

    `Vector2` is implicitly convertible to `std::complex<double>`.

??? func "`#!cpp Vector2::operator<<()`"

    `Vector2` can be serialized.

    ```cpp
    Vector2 v{1.2, 3.4};
    std::cout << v << std::endl;
    // prints something like: <1.2, 3.4>
    ```

### Arithmetic

Vector2 supports the element-wise addition, subraction, and scalar multiplication you would probably expect.

#### Rotations and complex multiplication

Our `Vector2` types further obey the multiplication and division rules of complex arithmetic, and thus can be used to represent rotations. For instance, a unit 2D vector representing a rotation can be used to rotate another vector like:
```cpp
Vector2 v = /* your vector */
Vector2 r = Vector2::fromAngle(PI/4); // rotation by 45 degrees
Vector2 vRot = r * v;
```
This is fundamentally no different from using 2x2 rotation matrices, but leads to much cleaner code (try using division to compute relative rotations!).


### Member operations

These methods _do not_ change the underlying `Vector2`, but return a new `Vector2`.
```cpp
Vector2 vec{1., 2.};
vec.rotate90();         // does nothing
vec = vec.rotate90();   // much better
```

??? func "`#!cpp Vector2 Vector2::normalize()`"

    Returns a unit-norm vector with the same direction. If the input is the zero vector, the result will contain NaNs.

??? func "`#!cpp Vector2 Vector2::normalizeCutoff(double mag = 0.)`"

    Returns a unit-norm vector with the same direction. If the input has magnitude less `<= mag`, the vector is unchanged.

??? func "`#!cpp Vector2 Vector2::unit()`"

    Alias for `normalize()`. 

??? func "`#!cpp Vector2 Vector2::rotate(double theta)`"

    Rotate the vector by angle $\theta$ in the counter-clockwise direction.

??? func "`#!cpp Vector2 Vector2::rotate90()`"

    Rotate the vector by $90^{\circ}$ in the counter-clockwise direction.

??? func "`#!cpp Vector2 Vector2::pow(double p)`"

    Raise the vector to a real power, in the sense of complex arithmetic. (see [std::pow](https://en.cppreference.com/w/cpp/numeric/complex/pow))

??? func "`#!cpp Vector2 Vector2::pow(Vector2 p)`"

    Raise the vector to a complex power, in the sense of complex arithmetic. (see [std::pow](https://en.cppreference.com/w/cpp/numeric/complex/pow))

??? func "`#!cpp Vector2 Vector2::conj()`"

    Transform the vector to its complex conjugate, negating the `y` component.

??? func "`#!cpp Vector2 Vector2::inv()`"

    Invert the vector, in the sense of complex arithmetic. Equivalent to `Vector2{1., 0.} / v`.


### Function operations

These operations do not change the vector on which they are called.

??? func "`#!cpp double norm(Vector2 v)`"

    Returns the magnitude of the vector.

    Also available as `v.norm()`.


??? func "`#!cpp double norm2(Vector2 v)`"

    Returns the squared magnitude of the vector.

    Also available as `v.norm()`.


??? func "`#!cpp Vector2 normalize(Vector2 v)`"

    Returns normalized copy of the vector.

??? func "`#!cpp Vector2 normalizeCutoff(Vector2 v, double mag = 0.)`"

    Returns a normalized copy of the vector. If the input has magnitude less `<= mag`, the vector is unchanged.

??? func "`#!cpp Vector2 unit(Vector2 v)`"

    Alias for `normalize(v)`.

??? func "`#!cpp double arg(Vector2 v)`"

    Returns the argument in the sense of complex arithmetic (i.e., the angle against the $x$-axis).

    Also available as `v.arg()`.


??? func "`#!cpp double dot(Vector2 u, Vector2 v)`"

    Returns the dot product between two vectors.


??? func "`#!cpp double cross(Vector2 u, Vector2 v)`"

    Returns the "cross" product between two vectors, that is `u.x * v.y - u.y * v.x`. Intuitively, the $z$-component of the 3D cross product of vectors in the plane.


??? func "`#!cpp Vector3 cross3(Vector2 u, Vector2 v)`"

    Returns the 3D cross product of vectors in the plane.


??? func "`#!cpp double angle(Vector2 u, Vector2 v)`"

    Returns the angle between two not-necessarily-unit vectors. Output in the range $[0, \pi]$.

??? func "`#!cpp Vector2 clamp(Vector2 val, Vector2 low, Vector2 high)`"

    Returns returns a a vector where each component has been clamped to be between the corresponding compnents of `low` and `high`.


??? func "`#!cpp Vector2 componentwiseMin(Vector2 u, Vector2 v)`"

    Returns a new vector, each component of which is the minimum of that component in `u` and `v`.


??? func "`#!cpp Vector2 componentwiseMax(Vector2 u, Vector2 v)`"

    Returns a new vector, each component of which is the maximum of that component in `u` and `v`.


### Properties

??? func "`#!cpp bool isfinite(Vector2 u)`"

    Returns true if both of the components of the vector are finite.

    Note: this function is intentionally not camel-cased out of solidarity with `std::isfinite()`.

    Also available as `u.isFinite()`.


??? func "`#!cpp bool isDefined(Vector2 u)`"

    Returns true if both of the components of the vector are not NaN.
    
    Also available as `u.isDefined()`.
