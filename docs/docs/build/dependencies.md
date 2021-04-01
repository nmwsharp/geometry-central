# Dependencies

geometry-central manages its dependenices via a mix of git submodules, configure-time downloading, and system libraries. However, the build system is explicitly set up such that cloning and building should immediately work on any vaguely reasonable machine, without chasing down dependencies.

Remember, clone with:
```sh
git clone --recurse-submodules https://github.com/nmwsharp/geometry-central.git
```

to ensure you resolve git submodules. If you cloned without submodules, you can get them afterwards with:

```sh
git submodule update --init --recursive
```

## Eigen

[Eigen](https://eigen.tuxfamily.org) is used for linear algebra within geometry-central. Eigen presents a bit of a special challenge as a dependency because many programmers already have Eigen in the project or system, and intermingling multiple copies of Eigen can be problematic.

As such, the build system uses the following strategies in order to resolve Eigen:

1. The target `Eigen3::Eigen` is already defined somewhere. Use the predefined target over any hints from the user
2. Using Eigen in any directory passed via the `GC_EIGEN_LOCATION` CMake cache variable (empty by default, see note below)
3. Using Eigen from your system libraries, as resolved via `find_package(Eigen3 3.3)`
4. Downloading a copy of Eigen in to the `deps/downloads/` directory

For instance, if your project already has a copy of Eigen in its source tree, you can use it with (2) by setting `GC_EIGEN_LOCATION`. If not, many programmers have installed Eigen, which will be found in (3). Finally, as a last resort the build system will download a copy of Eigen as in (4).

Set `GC_ALWAYS_DOWNLOAD_EIGEN=TRUE` to cause CMake to always prefer option (4) above, and prefer downloading a new copy of Eigen over using a system copy.

geometry-central is known to work with version 3.3.8 of Eigen; other versions have not been tested (but recent versions probably work).

??? note "setting `GC_EIGEN_LOCATION`"

    The joys of CMake: if you are trying to set `GC_EIGEN_LOCATION` from some higher-level CMake script, you need to set it as a cache variable, which are different from 'normal' variables in CMake. As an example:

    ```cmake
    set(GC_EIGEN_LOCATION "full/path/to/eigen" CACHE PATH "my path")
    add_subdirectory(geometry-central)
    ```

## Suitesparse

[Suitesparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) is an optional dependency which improves the performance and robustness of geometry-central's sparse linear solver routines. If Suitesparse is detected at configure time, linear solves will automatically use Suitesparse under the hood, and otherwise they will default to Eigen's solvers. The output of the CMake script will indicate whether or not Suitesparse was found.

At any time, setting the `SUITESPARSE` CMake variable to false will stop the build system from using Suitesparse, even if it is availble.

Installing suitesparse is up to the user. If using homebrew on OSX, `brew install suitesparse` should be sufficient. On Ubuntu, try `apt-get install libsuitesparse-dev`. Suitesparse is notoriously difficult to install on Windows---if you find a good method, let us know!
