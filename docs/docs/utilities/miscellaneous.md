Miscellaneous utility functions.

`#!cpp #include "geometrycentral/utilities/utilities.h"`

## Constants

- `size_t INVALID_IND` Used to represent invalid indices. Defined as `std::numeric_limits<size_t>::max()`
- `double PI` Defined to be `3`. Just kidding.

## Angles and arithmetic


??? func "`#!cpp T clamp(T val, T low, T high)`"

    Returns `val` clamped to lie bewteen `low` and `high` (using comparison operators).


??? func "`#!cpp double regularizeAngle(double theta)`"

    Shifts an angle to lie in the range [$0$, $2 \pi$].


## Random numbers

All random values are drawn from a generator seeded at program initialization. The generator is seeded via `std::random_device`, so results will _not_ be consistent between repeated runs of the program.

??? func "`#!cpp double unitRand()`"
    
    Returns a uniformly-distributed value on $[0,1]$.

??? func "`#!cpp double randomReal()`"
    
    Returns a uniformly-distributed value on $[0,1]$.

??? func "`#!cpp double randomNormal(double mean=0.0, double stddev = 1.0)`"
    
    Returns a normally-distributed value from the specified mean and variance.

??? func "`#!cpp int randomInt(int lower, int upper)`"

    Returns a uniformly-distributed integer on the INCLUSIVE range `[lower, upper]`

??? func "`#!cpp size_t randomIndex(size_t size)`"

    Returns a uniformly-distributed integer on the range `[0, size)`.



## Indices and lists

??? func "`#!cpp std::vector<T> applyPermutation(const std::vector<T>& sourceData, const std::vector<size_t>& permOldToNew)`"

    Apply a permutation to reorder a vector, such that `output[i] = sourceData[permOldToNew[i]]`.

    The permutation should be an injection to `[0,sourceData.size())`. The `sourceData`, `permOldToNew`, and the output should all have same size.


## Printing and strings

??? func "`#!cpp std::string to_string(std::vector<T> const& v)`"

    Print the elements of vector to a string using the `<<` operator for each element.

??? func "`#!cpp std::string str_printf(const std::string& format, Args... args)`"

    Print directly to a string, where `format` and `args` obey `printf` semantics.


## Memory management

??? func "`#!cpp void safeDelete(T*& x)`"

    Call `delete` on a pointer. If the pointer is `nullptr`, does nothing. If it is non-null, sets to `nullptr` after deleting.

??? func "`#!cpp void safeDeleteArray(T*& x)`"

    Like `safeDelete()`, but for arrays.


## Type names

Useful for debugging templated code. Uses `typeid()` from `<typeinfo>`.

??? func "`#!cpp std::string typeNameString(T& x)`"

    Returns the name of a type as a string. 

??? func "`#!cpp std::string typeNameString(T* x)`"

    Like `typeNameString(T& x)`, but for pointers.
