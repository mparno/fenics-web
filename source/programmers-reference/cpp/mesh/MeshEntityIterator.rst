.. Documentation for the header file dolfin/mesh/MeshEntityIterator.h

.. _programmers_reference_cpp_mesh_meshentityiterator:

MeshEntityIterator.h
====================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: MeshEntityIterator

    MeshEntityIterator provides a common iterator for mesh entities
    over meshes, boundaries and incidence relations. The basic use
    is illustrated below.
    
    The following example shows how to iterate over all mesh entities
    of a mesh of topological dimension dim:
    
        for (MeshEntityIterator e(mesh, dim); !e.end(); ++e)
        {
          e->foo();
        }
    
    The following example shows how to iterate over mesh entities of
    topological dimension dim connected (incident) to some mesh entity f:
    
        for (MeshEntityIterator e(f, dim); !e.end(); ++e)
        {
          e->foo();
        }
    
    In addition to the general iterator, a set of specific named iterators
    are provided for entities of type Vertex, Edge, Face, Facet and Cell.
    These iterators are defined along with their respective classes.

    .. cpp:function:: //    MeshEntityIterator(const MeshEntityIterator& entity) :  entity(entity.entity.mesh(), 0, 0), _pos(0)
                       //
    
        Copy constructor is private to disallow usage. If it were public (or not
        declared and thus a default version available) it would allow code like
        
        for (CellIterator c0(mesh); !c0.end(); ++c0)
          for (CellIterator c1(c0); !c1.end(); ++c1)
             ...
        
        c1 looks to be an iterator over the entities around c0 when it is in
        fact a copy of c0.

    .. cpp:function:: //dolfin/mesh/MeshEntityIterator.h:94: Warning|508| Declaration of 'operator ==' shadows declaration accessible via operator->(),
                                                                                                                                                     //Use const_cast to use operator* inside comparison, which automatically
                                                                                                                                                     //updates the entity index corresponding to pos *before* comparison (since
                                                                                                                                                     //update of entity delays until request for entity)
                                                                                                                                                     bool operator==(const MeshEntityIterator & it) const
    
        Comparison operator.
        @internal
        Uncommenting following  results into the warning message:

    .. cpp:function:: //std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: MeshEntity& operator*()
    
        Dereference operator

    .. cpp:function:: MeshEntity& operator[] (uint pos)
    
        Random access operator.

    .. cpp:function:: MeshEntity* operator->()
    
        Member access operator

    .. cpp:function:: MeshEntityIterator end_iterator()
    
        Provide a safeguard iterator pointing beyond the end of an iteration
        process, either iterating over the mesh /or incident entities. Added to
        be bit more like STL iteratoren, since many algorithms rely on a kind of
        beyond iterator.

    .. cpp:function:: MeshEntityIterator& operator++()
    
        Step to next mesh entity (prefix increment)

    .. cpp:function:: MeshEntityIterator& operator--()
    
        Step to the previous mesh entity (prefix decrease)

    .. cpp:function:: MeshEntityIterator() : _pos(0), pos_end(0), index(0)
    
        Default constructor

    .. cpp:function:: MeshEntityIterator(const Mesh& mesh, uint dim)
                                         : entity(mesh, dim, 0), _pos(0), pos_end(mesh.size(dim)), index(0)
    
        Create iterator for mesh entities over given topological dimension

    .. cpp:function:: MeshEntityIterator(const MeshEntity& entity, uint dim)
                                         : entity(entity.mesh(), dim, 0), _pos(0), index(0)
    
        Create iterator for entities of given dimension connected to given entity

    .. cpp:function:: MeshEntityIterator(const MeshEntityIterator& it) :  entity(it.entity),
                                         _pos(it._pos), pos_end(it.pos_end), index(it.index)
    
        Copy Constructor

    .. cpp:function:: bool end() const
    
        Check if iterator has reached the end

    .. cpp:function:: uint pos() const
    
        Return current position

    .. cpp:function:: virtual ~MeshEntityIterator()
    
        Destructor

    .. cpp:function:: void set_end()
    
        Set pos to end position. To create a kind of mesh.end() iterator.

