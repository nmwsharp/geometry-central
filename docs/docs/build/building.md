# Building

geometry-central uses CMake for to configure the build system. The basic workflow for downloading and compiling geometry-central via a terminal is:

```sh
git clone --recurse-submodules https://github.com/nmwsharp/geometry-central.git
cd geometry-central
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

However, since geometry-central is just a library, this does not build any executables, it merely compiles the library.

You can add geometry-central to an existing project's `CMakeLists.txt` like

```cmake
add_subdirectory("path/to/geometry-central") # wherever you put it
target_link_libraries(your-project-target geometry-central)
```

#### Example

For a simple example project using geometry-central (with [Polyscope](https://polyscope.run) for visualization), see [gc-polyscope-project-template](https://github.com/nmwsharp/gc-polyscope-project-template). This is a good starting point for new projects using geometry-central.

### On Windows

When using Visual Studio on Windows, CMake can be used (either via the terminal or gui) to generate Visual Studio project and solution files. The project has been verified to compile out of the box with Visual Studio 2017 & 2019 (older versions not tested).

## Compile flags & options

The library includes a few optional safety checks which are performed at runtime, even in release mode. Such checks are generally very cheap yet quite useful. Nonetheless, adding the `NGC_SAFETY_CHECKS` define will disable all optional safety checks, for a very small increase in performance.

## Bazel

geometry-central is now available as a Bazel module: https://registry.bazel.build/modules/geometry-central

To use it, simply add this line to your MODULE.bazel:

`bazel_dep(name = "geometry-central", version = "1.0.0")`

and then specify the dependencies in your build rules. There are four public targets, which correspond to the groupings in the source code: 
 - `@geometry-central//:numerical`
 - `@geometry-central//:surface`
 - `@geometry-central//:pointcloud`
 - `@geometry-central//:utilities`

For example:

```python
cc_binary(
    name = "main",
    srcs = ["main.cpp"],
    data = ["fox.ply"],
    deps = [
        "@geometry-central//:surface",
        "@geometry-central//:utilities",
    ],
)
```