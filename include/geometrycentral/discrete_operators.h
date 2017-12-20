#pragma once

#include "geometry.h"

#include "Eigen/SparseCore"

namespace geometrycentral {

/*
    === Standard operators from DDG, using circumcentric dual, defined over the
   entire mesh ===

    All use the standard indexing convention according to
   mesh->getVertexIndices() (etc).
    Edge orientation convention is according to the direction of edge.halfedge.

    Note that discrete_operators.cpp contains explicit instantiations of these
   for T={double,Complex},
    and D=Euclidean. We should (TODO) probably move that to a .ipp file so they
   can be instantiated for
    user-specified types, then use macro-magic to explicitly instantiate them
   for the popular types.
    Right now, they effectively same real-valued matrix is being built in all
   cases, its just the type
    that differs.
*/

// TODO Implement a combined operator that builds all (many?) of these at once,
// saving recomputation

// Hodge star on 0-forms. Returns a (nVerts, nVerts) matrix.
Eigen::SparseMatrix<double> buildHodge0(Geometry<Euclidean>* geometry);

// Hodge star on 1-forms. Returns a (nEdges, nEdges) matrix.
Eigen::SparseMatrix<double> buildHodge1(Geometry<Euclidean>* geometry);

// Hodge star on 2-forms. Returns a (nFaces, nFaces) matrix.
Eigen::SparseMatrix<double> buildHodge2(Geometry<Euclidean>* geometry);

// Derivative on 0-forms. Returns a (nEdges, nVerts) matrix
Eigen::SparseMatrix<double> buildDerivative0(HalfedgeMesh* mesh);

// Derivative on 1-forms. Returns a (nFaces, nEdges) matrix
Eigen::SparseMatrix<double> buildDerivative1(HalfedgeMesh* mesh);


}  // namespace geometrycentral