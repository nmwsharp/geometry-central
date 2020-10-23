
// === Helpers which will allow us abstract over types
// The corresponding template declarations are given in halfedge_element_types.h

namespace geometrycentral {

// clang-format off

template<> inline size_t nElements<surface::Vertex       >(surface::SurfaceMesh* mesh)   { return mesh->nVertices();   }
template<> inline size_t nElements<surface::Face         >(surface::SurfaceMesh* mesh)   { return mesh->nFaces();      }
template<> inline size_t nElements<surface::Edge         >(surface::SurfaceMesh* mesh)   { return mesh->nEdges();      }
template<> inline size_t nElements<surface::Halfedge     >(surface::SurfaceMesh* mesh)   { return mesh->nHalfedges();  }
template<> inline size_t nElements<surface::Corner       >(surface::SurfaceMesh* mesh)   { return mesh->nCorners();    }
template<> inline size_t nElements<surface::BoundaryLoop >(surface::SurfaceMesh* mesh)   { return mesh->nBoundaryLoops();    }

template<> inline size_t elementCapacity<surface::Vertex      >(surface::SurfaceMesh* mesh)   { return mesh->nVerticesCapacity();   }
template<> inline size_t elementCapacity<surface::Face        >(surface::SurfaceMesh* mesh)   { return mesh->nFacesCapacity(); }
template<> inline size_t elementCapacity<surface::Edge        >(surface::SurfaceMesh* mesh)   { return mesh->nEdgesCapacity();      }
template<> inline size_t elementCapacity<surface::Halfedge    >(surface::SurfaceMesh* mesh)   { return mesh->nHalfedgesCapacity();}
template<> inline size_t elementCapacity<surface::Corner      >(surface::SurfaceMesh* mesh)   { return mesh->nHalfedgesCapacity();    }
template<> inline size_t elementCapacity<surface::BoundaryLoop>(surface::SurfaceMesh* mesh)   { return mesh->nBoundaryLoopsCapacity();    }

template<> inline size_t dataIndexOfElement<surface::Vertex          >(surface::SurfaceMesh* mesh, surface::Vertex e)         { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<surface::Face            >(surface::SurfaceMesh* mesh, surface::Face e)           { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<surface::Edge            >(surface::SurfaceMesh* mesh, surface::Edge e)           { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<surface::Halfedge        >(surface::SurfaceMesh* mesh, surface::Halfedge e)       { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<surface::Corner          >(surface::SurfaceMesh* mesh, surface::Corner e)         { return e.getIndex(); }
template<> inline size_t dataIndexOfElement<surface::BoundaryLoop    >(surface::SurfaceMesh* mesh, surface::BoundaryLoop e)   { return e.getIndex(); }

template<> inline surface::VertexSet         iterateElements<surface::Vertex      >(surface::SurfaceMesh* mesh)   { return mesh->vertices();      }
template<> inline surface::HalfedgeSet       iterateElements<surface::Halfedge    >(surface::SurfaceMesh* mesh)   { return mesh->halfedges();     }
template<> inline surface::CornerSet         iterateElements<surface::Corner      >(surface::SurfaceMesh* mesh)   { return mesh->corners();       }
template<> inline surface::EdgeSet           iterateElements<surface::Edge        >(surface::SurfaceMesh* mesh)   { return mesh->edges();         }
template<> inline surface::FaceSet           iterateElements<surface::Face        >(surface::SurfaceMesh* mesh)   { return mesh->faces();         }
template<> inline surface::BoundaryLoopSet   iterateElements<surface::BoundaryLoop>(surface::SurfaceMesh* mesh)   { return mesh->boundaryLoops(); }

template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Vertex      >(surface::SurfaceMesh* mesh)   { return mesh->vertexExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Halfedge    >(surface::SurfaceMesh* mesh)   { return mesh->halfedgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Corner      >(surface::SurfaceMesh* mesh)   { return mesh->halfedgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Edge        >(surface::SurfaceMesh* mesh)   { return mesh->edgeExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::Face        >(surface::SurfaceMesh* mesh)   { return mesh->faceExpandCallbackList;   }
template<> inline std::list<std::function<void(size_t)>>& getExpandCallbackList<surface::BoundaryLoop>(surface::SurfaceMesh* mesh)   { return mesh->faceExpandCallbackList;   }

template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Vertex       >(surface::SurfaceMesh* mesh)   { return mesh->vertexPermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Halfedge     >(surface::SurfaceMesh* mesh)   { return mesh->halfedgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Corner       >(surface::SurfaceMesh* mesh)   { return mesh->halfedgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Edge         >(surface::SurfaceMesh* mesh)   { return mesh->edgePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::Face         >(surface::SurfaceMesh* mesh)   { return mesh->facePermuteCallbackList;   }
template<> inline std::list<std::function<void(const std::vector<size_t>&)>>& getPermuteCallbackList<surface::BoundaryLoop >(surface::SurfaceMesh* mesh)   { return mesh->boundaryLoopPermuteCallbackList;   }

template<> inline std::string typeShortName<surface::Vertex       >()            { return "v";    }
template<> inline std::string typeShortName<surface::Halfedge     >()            { return "he";   }
template<> inline std::string typeShortName<surface::Corner       >()            { return "c";    }
template<> inline std::string typeShortName<surface::Edge         >()            { return "e";    }
template<> inline std::string typeShortName<surface::Face         >()            { return "f";    }
template<> inline std::string typeShortName<surface::BoundaryLoop >()            { return "bl";   }

// clang-format on

} // namespace geometrycentral
