
// === Helpers which will allow us abstract over types
// The corresponding template declarations are given in halfedge_element_types.h

namespace geometrycentral {
namespace surface {

// clang-format off

template<> inline size_t nElements<Vertex       >(HalfedgeMesh* mesh)   { return mesh->nVertices();   }
template<> inline size_t nElements<Face         >(HalfedgeMesh* mesh)   { return mesh->nFaces();      }
template<> inline size_t nElements<Edge         >(HalfedgeMesh* mesh)   { return mesh->nEdges();      }
template<> inline size_t nElements<Halfedge     >(HalfedgeMesh* mesh)   { return mesh->nHalfedges();  }
template<> inline size_t nElements<Corner       >(HalfedgeMesh* mesh)   { return mesh->nCorners();    }
template<> inline size_t nElements<BoundaryLoop >(HalfedgeMesh* mesh)   { return mesh->nBoundaryLoops();    }

template<> inline size_t elementCapacity<Vertex      >(HalfedgeMesh* mesh)   { return mesh->nVerticesCapacity();   }
template<> inline size_t elementCapacity<Face        >(HalfedgeMesh* mesh)   { return mesh->nFacesCapacity() + mesh->nBoundaryLoops(); }
template<> inline size_t elementCapacity<Edge        >(HalfedgeMesh* mesh)   { return mesh->nEdgesCapacity();      }
template<> inline size_t elementCapacity<Halfedge    >(HalfedgeMesh* mesh)   { return mesh->nHalfedgesCapacity();}
template<> inline size_t elementCapacity<Corner      >(HalfedgeMesh* mesh)   { return mesh->nHalfedgesCapacity();    }
template<> inline size_t elementCapacity<BoundaryLoop>(HalfedgeMesh* mesh)   { return mesh->nBoundaryLoopsCapacity();    }

template<> inline size_t dataIndexOfElement<Vertex          >(HalfedgeMesh* mesh, Vertex e)         { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<Face            >(HalfedgeMesh* mesh, Face e)           { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<Edge            >(HalfedgeMesh* mesh, Edge e)           { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<Halfedge        >(HalfedgeMesh* mesh, Halfedge e)       { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<Corner          >(HalfedgeMesh* mesh, Corner e)         { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<BoundaryLoop    >(HalfedgeMesh* mesh, BoundaryLoop e)   { return e.getIndex(); }

template <> struct ElementSetType<Vertex        >   { typedef VertexSet       type; };
template <> struct ElementSetType<Face          >   { typedef FaceSet         type; };
template <> struct ElementSetType<Edge          >   { typedef EdgeSet         type; };
template <> struct ElementSetType<Halfedge      >   { typedef HalfedgeSet     type; };
template <> struct ElementSetType<Corner        >   { typedef CornerSet       type; };
template <> struct ElementSetType<BoundaryLoop  >   { typedef BoundaryLoopSet type; };

template<> inline VertexSet         iterateElements<Vertex      >(HalfedgeMesh* mesh)   { return mesh->vertices();      }
template<> inline HalfedgeSet       iterateElements<Halfedge    >(HalfedgeMesh* mesh)   { return mesh->halfedges();     }
template<> inline CornerSet         iterateElements<Corner      >(HalfedgeMesh* mesh)   { return mesh->corners();       }
template<> inline EdgeSet           iterateElements<Edge        >(HalfedgeMesh* mesh)   { return mesh->edges();         }
template<> inline FaceSet           iterateElements<Face        >(HalfedgeMesh* mesh)   { return mesh->faces();         }
template<> inline BoundaryLoopSet   iterateElements<BoundaryLoop>(HalfedgeMesh* mesh)   { return mesh->boundaryLoops(); }

template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<Vertex      >(HalfedgeMesh* mesh)   { return mesh->vertexExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<Halfedge    >(HalfedgeMesh* mesh)   { return mesh->halfedgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<Corner      >(HalfedgeMesh* mesh)   { return mesh->halfedgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<Edge        >(HalfedgeMesh* mesh)   { return mesh->edgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<Face        >(HalfedgeMesh* mesh)   { return mesh->faceExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<BoundaryLoop>(HalfedgeMesh* mesh)   { return mesh->faceExpandCallbackList;   }

template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<Vertex       >(HalfedgeMesh* mesh)   { return mesh->vertexPermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<Halfedge     >(HalfedgeMesh* mesh)   { return mesh->halfedgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<Corner       >(HalfedgeMesh* mesh)   { return mesh->halfedgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<Edge         >(HalfedgeMesh* mesh)   { return mesh->edgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<Face         >(HalfedgeMesh* mesh)   { return mesh->facePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<BoundaryLoop >(HalfedgeMesh* mesh)   { return mesh->facePermuteCallbackList;   }

template<> inline std::string typeShortName<Vertex       >()            { return "v";    }
template<> inline std::string typeShortName<Halfedge     >()            { return "he";   }
template<> inline std::string typeShortName<Corner       >()            { return "c";    }
template<> inline std::string typeShortName<Edge         >()            { return "e";    }
template<> inline std::string typeShortName<Face         >()            { return "f";    }
template<> inline std::string typeShortName<BoundaryLoop >()            { return "bl";   }

// clang-format on

} // namespace surface
} // namespace geometrycentral
