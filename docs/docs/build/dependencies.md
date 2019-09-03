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

1. Using Eigen in any directory passed via the `GC_EIGEN_LOCATION` CMake variable (empty by default)
2. Using Eigen from your system libraries, as resolved via `find_package(Eigen3 3.3)`
3. Downloading a copy of Eigen in to the `deps/downloads/` directory

For instance, if your project already has a copy of Eigen in its source tree, you can use with (1) by setting `GC_EIGEN_LOCATION`. If not, many programmers have installed Eigen, which will be found in (2). Finally, as a last resort the build system will download a copy of Eigen as in (3).

geometry-central is known to work with version 3.3 of Eigen; other versions have not been tested (but recent versions probably work).

Note: once upon a time, Eigen was a submodule of geometry-central. If updating from an old version, you may need to manually delete `deps/eigen-git-mirror`.

## Suitesparse

geometry-central's linear solvers will automatically use [Suitesparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) routines under the hood if detected at configure time. The output of the CMake script will indicate whether or not suitesparse was found.

At any time, setting the `SUITESPARSE` CMake variable to false will stop the build system from using Suitesparse, even if it is availble.

Installing suitesparse is up to the user. If using homebrew on OSX, `brew install suitesparse` should be sufficient. On Ubuntu, try `apt-get install libsuitesparse-dev`. Suitesparse is notoriously difficult to install on Windows---if you find a good method, let us know!
