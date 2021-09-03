`#include "geometrycentral/surface/transfer_functions.h"`

[Intrinsic Triangulations](/surface/intrinsic_triangulations/basics) provide high quality _function spaces_, even on near-degenerate geometry, which dramatically improves the accuracy of PDE-based algorithms such as the [vector heat method](/surface/algorithms/vector_heat_method). However, these accurate solutions live on the intrinsic triangulation, and cannot be represented exactly as piecewise-linear functions on the original mesh. Geometry-central provides several ways of transferring functions between different triangulations.

### Pointwise Transfer
The simplest way to transfer functions between two triangulations of the same surface is to simply evaluate the function at the vertices of the second triangulation and then interpolate these values piecewise-linearly.

### L2-Optimal Transfer
Rather than sampling values, one can instead look for the piecewise linear function on the second triangulation which is closest in the $L^2$-sense to the target function on the first triangulation. This can be computed exactly via a linear least squares system. See [Integer Coordinates for Intrinsic Geometry Processing](https://arxiv.org/pdf/2106.00220.pdf) for more details.

## Single Transfer Functions
A one-off utility function is provided which transfers functions between different triangulations of the same surface.
Repeated solves should use the stateful version below.

Example
```cpp
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/transfer_functions.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Create an intrinsic triangulation
SignpostIntrinsicTriangulation intTri(*mesh, *geometry);

// Change the intrinsic triangulation
intTri.delaunayRefine();

// Compute something useful on the intrinsic triangulation
VertexData<double> f_intrinsic = /* some function */

// Transfer function back to extrinsic mesh
VertexData<double> f_extrinsic = transferBtoA(intTri, f_intrinsic, TransferMethod::L2);
```

??? func "`#!cpp VertexData<double> transferAtoB(IntrinsicTriangulation& intTri, const VertexData<double>& valuesOnA, TransferMethod method)`"
    Transfers a scalar function from `intTri.inputMesh` to `intTri.intrinsicMesh`
    
    - `valuesOnA` : the data on `intTri.inputMesh` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.

??? func "`#!cpp VertexData<double> transferBtoA(IntrinsicTriangulation& intTri, const VertexData<double>& valuesOnB, TransferMethod method)`"
    Transfers a scalar function from `intTri.intrinsicMesh` to `intTri.inputMesh`
    
    - `valuesOnB` : the data on `intTri.intrinsicMesh` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.
    

The `method` argument is an enum:
```cpp
enum class TransferMethod { Pointwise = 0, L2 };
```
    
    
## Repeated Transfer Functions

The stateful class `AttributeTransfer` transfers functions between different triangulations.

```cpp
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/signpost_intrinsic_triangulation.h"
#include "geometrycentral/surface/transfer_functions.h"

// Load a mesh
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::tie(mesh, geometry) = loadMesh(filename);

// Create an intrinsic triangulation
SignpostIntrinsicTriangulation intTri(*mesh, *geometry);

// Change the intrinsic triangulation
intTri.delaunayRefine();

// Create the AttributeTransfer object
AttributeTransfer transfer(intTri);

// Compute several functions on the intrinsic triangulation
std::vector<VertexData<double>> intrinsicFunctions = /* some functions */

// Transfer functions back to extrinsic mesh
for (VertexData<double> f_intrinsic : intrinsicFunctions) {
    VertexData<double> f_extrinsic = transfer.transferBtoA(f_intrinsic, TransferMethod::L2);
    /* do something with f_extrinsic */
}

```

### Constructors

??? func "`#!cpp AttributeTransfer::AttributeTransfer(CommonSubdivision& cs, VertexPositionGoemetry& geomA)`"

    Create a new solver for attribute transfer. Precomputation is performed lazily as needed.
    
    - `cs` is the common subdivision of the two triangulations between which the solver will transfer data.
    
    - `geomA` is the geometry of `meshA`.
    
??? func "`#!cpp AttributeTransfer::AttributeTransfer(CommonSubdivision& cs, IntrinsicTriangulation& intTri)`"

    Create a new solver for attribute transfer. Precomputation is performed lazily as needed.
    
    - `cs` is the common subdivision of the two triangulations between which the solver will transfer data.
    
    - `intTri` is the `IntrinsicTriangulation` object representing the two triangulations.
    
### Methods
??? func "`#!cpp VertexData<double> AttributeTransfer::transferAtoB(const VertexData<double>& valuesOnA, TransferMethod method)`"
    Transfers a scalar function from `meshA` to `meshB`
    
    - `valuesOnA` : the data on `meshA` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.

??? func "`#!cpp VertexData<double> AttributeTransfer::transferBtoA(const VertexData<double>& valuesOnB, TransferMethod method)`"
    Transfers a scalar function from `meshB` to `meshA`
    
    - `valuesOnB` : the data on `meshB` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.
