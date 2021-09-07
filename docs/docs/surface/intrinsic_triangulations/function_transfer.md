[Intrinsic triangulations](/surface/intrinsic_triangulations/basics) provide high quality _function spaces_, even on near-degenerate geometry, which dramatically improves the accuracy of PDE-based algorithms (e.g. the heat-based methods in geometry-central). However, these accurate solutions live on the intrinsic triangulation, and cannot be represented exactly as piecewise-linear functions on the original mesh. 

Given an input mesh of some surface, and an [intrinsic triangulations](/surface/intrinsic_triangulations/basics) of that surface (which may have additional vertices inserted), these functions allow you to transfer values between the two domains via several different strategies. 

## Pointwise Transfer at Vertices

The simplest method to transfer functions between two triangulations of the same surface is to simply copy values at vertices.  The following methods, from the `IntrinsicTriangulation` class, implement this simple transfer scheme in both the forward and reverse direction.

??? func "`#!cpp VertexData<T> IntrinsicTriangulation::sampleFromInput(const VertexData<T>& dataOnInput)`"
  
    Given data defined on the vertices of the input triangulation, samples it to the vertices of the intrinsic triangulation.

    If the intrinsic triangulation contains new, inserted vertices which are not in common with the input triangulation, they sample a linearly-interpolated value from the input function.
   

??? func "`#!cpp VertexData<T> IntrinsicTriangulation::restrictToInput(const VertexData<T>& dataOnIntrinsic)`"
  
    Given data defined on the vertices of the intrinsic triangulation, restrict it to the vertices of the input triangulation (that is, simply copy values from the shared vertices).
    
    If the intrinsic triangulation contains new, inserted vertices which are not in common with the input triangulation, their values are ignored for the purposes of restriction.



## Transfer to the Common Subdivision

Another possibility is to transfer data to the [common subdivison](/surface/intrinsic_triangulations/common_subdivision), a special triangulation which is a superset of both the input and the intrinsic triangulation. This transfer is easy to define, because the common subdivision is precisely the triangulation whose linear bases can exactly represent functions from either the input or intrinsic triangulation.

Of course, unlike the other methods described on this page, this strategy _does not_ transfer functions directly between an input triangulation and an intrinsic triangulation, but rather transfers values from either the input or intrinsic triangulation _to_ the common subdivision. 

One common use for transferring values to the common subdivision is visualization, because the common subdivision is naturally embedded in space as an ordinary mesh.

??? func "`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossA(const VertexData<T>& dataA)`<br/>`#!cpp VertexData<T> CommonSubdivision::interpolateAcrossB(const VertexData<T>& dataB)`"

    Linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision.


??? func "`#!cpp FaceData<T> CommonSubdivision::copyFromA(const FaceData<T>& dataA)`<br/>`#!cpp FaceData<T> CommonSubdivision::copyFromB(const FaceData<T>& dataB)`"

    Copy data at faces from one of the meshes to the common subdivision. Each face of the common subdivision gets the value from the face which contains it. The return value is defined per-face of the common subdivision mesh.

??? func "`#!cpp SparseMatrix<double> CommonSubdivison::interpolationMatrixA()`<br/>`#!cpp SparseMatrix<double> CommonSubdivision::interpolationMatrixB()`"

    Yields a `|V| x |V_A|` matrix (resp. `|V| x |V_B|`) which linearly interpolates data on `meshA` (resp. `meshB`) to the common subdivision. Here `|V|` denotes the number of vertices in the common subdivision.


## L2-Optimal Transfer

`#include "geometrycentral/surface/transfer_functions.h"`

When transferring data directly between the input triangulation and intrinsic triangulation, simple sampling values is naive, and has no reason to be the "best" approach.  Instead, one can directly compute the function on the other surface, in the sense of $L_2$-distance between functions---see [Integer Coordinates for Intrinsic Geometry Processing](https://arxiv.org/pdf/2106.00220.pdf) for details.

Computationally, this amounts to solving a sparse linear least-squares problem defined via the common subdivision. Though somewhat more expensive, this approach can greatly improve accuracy.

### Single Transfer Functions

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
    

### Transfer Method

The `method` argument is an enum:

```cpp
enum class TransferMethod { Pointwise = 0, L2 };
```

for completeness, the`TransferMethod::Pointwise` option implements the [pointwise transfer](/surface/intrinsic_triangulations/function_transfer/#pointwise-transfer-at-vertices) described above, while `TransferMethod::L2` is the optimal scheme described here.

    
    
### Repeated Transfer Functions

If many functions are to be transferred, pre-factoring the linear least-squares problem can greatly improve performance.  The stateful class `AttributeTransfer` facilitates this precomputation.

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

#### Constructors

??? func "`#!cpp AttributeTransfer::AttributeTransfer(CommonSubdivision& cs, VertexPositionGoemetry& geomA)`"

    Create a new solver for attribute transfer. Precomputation is performed lazily as needed.
    
    - `cs` is the common subdivision of the two triangulations between which the solver will transfer data.
    
    - `geomA` is the geometry of `meshA`.
    
??? func "`#!cpp AttributeTransfer::AttributeTransfer(CommonSubdivision& cs, IntrinsicTriangulation& intTri)`"

    Create a new solver for attribute transfer. Precomputation is performed lazily as needed.
    
    - `cs` is the common subdivision of the two triangulations between which the solver will transfer data.
    
    - `intTri` is the `IntrinsicTriangulation` object representing the two triangulations.
    
#### Methods

??? func "`#!cpp VertexData<double> AttributeTransfer::transferAtoB(const VertexData<double>& valuesOnA, TransferMethod method)`"
    Transfers a scalar function from `meshA` to `meshB`
    
    - `valuesOnA` : the data on `meshA` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.

??? func "`#!cpp VertexData<double> AttributeTransfer::transferBtoA(const VertexData<double>& valuesOnB, TransferMethod method)`"
    Transfers a scalar function from `meshB` to `meshA`
    
    - `valuesOnB` : the data on `meshB` to be transferred.
    
    - `method` : either `TransferMethod::Pointwise` for pointwise transfer of `TransferMethod::L2` for $L^2$-optimal transfer.
