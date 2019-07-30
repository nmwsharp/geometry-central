# Linear algebra utilities

### Construct and convert

??? func "`#!cpp SparseMatrix<T> identityMatrix(size_t N)`"

    Construct and `N x N` identity matrix of the requested type.


??? func "`#!cpp void shiftDiagonal(SparseMatrix<T>& m, T shiftAmount = 1e-4)`"

    Shift the diagonal of matrix, by adding `A + shiftDiagonal * identityMatrix()`.


??? func "`#!cpp SparseMatrix<T> verticalStack(const std::vector<SparseMatrix<T>>& mats)`"

    Vertically stack sparse matrices like

    $$
    A,B,C \to
    \begin{bmatrix}
    A \\
    B \\
    C
    \end{bmatrix}
    $$

    all matrices must have the same number of columns.

    Example:
    ```cpp
    SparseMatrix<double> matA = /* 35 x N */
    SparseMatrix<double> matB = /* 10 x N */
    SparseMatrix<double> matC = /* 2 x N */
    SparseMatrix<double> stacked = verticalStack<double>({matA, matB, matC});
    ```

??? func "`#!cpp SparseMatrix<T> horizontalStack(const std::vector<SparseMatrix<T>>& mats)`"

    Vertically stack sparse matrices like

    $$
    A,B,C \to
    \begin{bmatrix}
    A & B & C
    \end{bmatrix}
    $$

    all matrices must have the same number of rows..

    Example:
    ```cpp
    SparseMatrix<double> matA = /* N x 35 */
    SparseMatrix<double> matB = /* N x 10 */
    SparseMatrix<double> matC = /* N x 2 */
    SparseMatrix<double> stacked = horizontalStack<double>({matA, matB, matC});
    ```

??? func "`#!cpp SparseMatrix<double> complexToReal(const SparseMatrix<std::complex<double>>& m)`"

    Convert an `N x M` complex matrix to a `2N x 2M` real matrix, expanding each complex component in to a `2 x 2` block to evaluate the complex product.


??? func "`#!cpp Vector<double> complexToReal(const Vector<std::complex<double>>& v)`"

    Convert an length `N` complex vector to a length `2N` real vector, expanding each complex component in to consecutive real and imaginary components.


### Validate matrix properties

??? func "`#!cpp void checkFinite(const Eigen::Matrix<>& m)`"

    Verify that all entries in an matrix are finite, throwing if not. Defined for all Eigen matrix, vector, and sparse matrix types.


??? func "`#!cpp void checkSymmetric(const Eigen::SparseMatrix<>& m, double absoluteEPS=-1.)`"

    Verify that a matrix is symmetric, throwing if not. Defined for all Eigen sparse matrix types.

    `absoluteEPS` is an epsilon to use for the element-wise comparison test. If the default value of `-1` is given, a reasonable epsilon is automatically computed from the matrix entries.

??? func "`#!cpp void checkHermitian(const Eigen::SparseMatrix<>& m, double absoluteEPS=-1.)`"

    Verify that a matrix is Hermitian, throwing if not. Defined for all Eigen sparse matrix types.

    For real matrices, identical to check symmetric.

    `absoluteEPS` is an epsilon to use for the element-wise comparison test. If the default value of `-1` is given, a reasonable epsilon is automatically computed from the matrix entries.


### Block decomposition

These routines assist with decomposing a square matrix in to interleaved submatrix blocks, where the blocks might not necessarily be contiguous. One common usage is extracting boundary components of a finite element matrix to apply boundary conditions, as in the example below.

Example usage:
```cpp

// Hypothetical input data
SparseMatrix<double> mat = /* your square matrix */;
size_t N = mat.rows();
size_t NBoundary = /* ... */;
Vector<double> rhsVals = Vector<double>::Zero(N);        // rhs for the system
Vector<double> bcVals = Vector<double>::Ones(NBoundary); // boundary values at
                                                         // some nodes

// Build the membership vector, which indicates which entries should be separated
// in to set "A" (others are in "B")
Vector<bool> setAMembership(N);
for(size_t i = 0; i < N; i++) {
  if(/* element i is boundary */) {
    setAMembership(i) = true;
  } else {
    setAMembership(i) = false;
  }
}

// Construct the decomposition 
BlockDecompositionResult<double> decomp = 
  blockDecomposeSquare(mat, setAMembership, true);

// The four sub-blocks of the matrix are now in
// decomp.AA, decomp.AB, decomp.BA, decomp.BB

// Split up the rhs vector
Vector<double> rhsValsA, rhsValsB;
decomposeVector(decomp, rhsVals, rhsValsA, rhsValsB);

// Solve problem
Vector<double> combinedRHS = rhsValsA - decomp.AB * bcVals;
Vector<double> Aresult = solve(decomp.AA, combinedRHS);

// Combine the two boundary conditions and interior solution to a full vector
Vector<double> result = reassembleVector(decomp, Aresult, bcVals);

```

??? func "`#!cpp BlockDecompositionResult<T> blockDecomposeSquare(const SparseMatrix<T>& m, const Vector<bool>& Aset, bool buildBuildBside = true)`"

    Build a block decomposition of a matrix.

??? func "`#!cpp void decomposeVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vec, Vector<T>& vecAOut, Vector<T>& vecBOut)`"

    Use an existing block decomposition to partition a vector.

??? func "`#!cpp Vector<T> reassembleVector(BlockDecompositionResult<T>& decomp, const Vector<T>& vecA, const Vector<T>& vecB)`"
    
    Use an existing block decomposition to build a vector from partitioned pieces.
