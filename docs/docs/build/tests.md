## Building and running tests

Compile and run the tests with:
```sh
cd test
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j12 && ./bin/geometry-central-test
```

## Tests organization

All tests are stored in the `test/` subdirectory, which is not touched by the usual build system.

### googletest

We use [googletest](https://github.com/google/googletest) as a testing framework. Most users do not need the tests, so rather than packing it as a git submodule which would unavoidably be cloned, the build system downloads a binary of the googletest when the tests are built---a consequence is that you must have a network connection to build tests for the first time (TODO: enable using a system install of googletest).

###  Assets

The `tests/assets/` directory contains a handful of input files for various tests. The absolute paths to these files are baked in to the test executable by the build system, so moving this directory after compiling tests may cause problems. The disk footprint of assets should be kept as small as possible since they are stored in the library repository.
