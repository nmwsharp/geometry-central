
inline std::array<Vector2, 4>
NormalCoordinateIntrinsicTriangulation::layoutDiamond(Halfedge iHe) const {

    // Conventions:
    //  - iHe points from vertex 2 to vertex 0, other vertices are numbered ccw
    //  - iHe is incident on face A, other is face B
    //  - halfedges within face are numbered CCW as A0, A1, A2 (etc),
    //    starting with iHe and twin(iHe)
    //  - When we lay out the triangle, p3 is at the origin and
    //    edge 3-0 is along the X-axis
    //  - flips is always ccw, so iHe points from vertex 3 --> 1 after

    // Gather index values
    Halfedge iHeA0 = iHe;
    Halfedge iHeA1 = iHeA0.next();
    Halfedge iHeA2 = iHeA1.next();
    Halfedge iHeB0 = iHe.twin();
    Halfedge iHeB1 = iHeB0.next();
    Halfedge iHeB2 = iHeB1.next();

    // Gather length values
    double l01 = intrinsicEdgeLengths[iHeA1.edge()];
    double l12 = intrinsicEdgeLengths[iHeA2.edge()];
    double l23 = intrinsicEdgeLengths[iHeB1.edge()];
    double l30 = intrinsicEdgeLengths[iHeB2.edge()];
    double l02 = intrinsicEdgeLengths[iHeA0.edge()];

    // Lay out the vertices of the diamond
    Vector2 p3{0., 0.};
    Vector2 p0{l30, 0.};
    Vector2 p2 = layoutTriangleVertex(
        p3, p0, l02, l23); // involves more arithmetic than strictly necessary
    Vector2 p1 = layoutTriangleVertex(p2, p0, l01, l12);

    return {p0, p1, p2, p3};
}
