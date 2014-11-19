
.. Documentation for the header file dolfin/mesh/MeshEntityIterator.h

.. _programmers_reference_cpp_mesh_meshentityiterator:

MeshEntityIterator.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEntityIterator

    MeshEntityIterator provides a common iterator for mesh entities
    over meshes, boundaries and incidence relations. The basic use
    is illustrated below.
    
    *Example*
    
        The following example shows how to iterate over all mesh entities
        of a mesh of topological dimension dim:
    
        .. code-block:: c++
    
            for (MeshEntityIterator e(mesh, dim); !e.end(); ++e)
            {
              e->foo();
            }
    
        The following example shows how to iterate over mesh entities of
        topological dimension dim connected (incident) to some mesh entity f:
    
        .. code-block:: c++
    
            for (MeshEntityIterator e(f, dim); !e.end(); ++e)
            {
              e->foo();
            }
    
    In addition to the general iterator, a set of specific named iterators
    are provided for entities of type :cpp:class:`Vertex`, :cpp:class:`Edge`, :cpp:class:`Face`, :cpp:class:`Facet`
    and :cpp:class:`Cell`. These iterators are defined along with their respective
    classes.


    .. cpp:function:: MeshEntityIterator()
    
        Default constructor


    .. cpp:function:: MeshEntityIterator(const Mesh& mesh, std::size_t dim)
    
        Create iterator for mesh entities over given topological dimension


    .. cpp:function:: MeshEntityIterator(const MeshEntity& entity, std::size_t dim)
    
        Create iterator for entities of given dimension connected to given entity


    .. cpp:function:: MeshEntityIterator(const MeshEntityIterator& it)
    
        Copy constructor


    .. cpp:function:: MeshEntityIterator& operator++()
    
        Step to next mesh entity (prefix increment)


    .. cpp:function:: MeshEntityIterator& operator--()
    
        Step to the previous mesh entity (prefix decrease)


    .. cpp:function:: std::size_t pos() const
    
        Return current position


    .. cpp:function:: bool operator==(const MeshEntityIterator& it) const
    
        Comparison operator


    .. cpp:function:: bool operator!=(const MeshEntityIterator & it) const
    
        Comparison operator


    .. cpp:function:: MeshEntity& operator*()
    
        Dereference operator


    .. cpp:function:: MeshEntity* operator->()
    
        Member access operator


    .. cpp:function:: MeshEntity& operator[] (std::size_t pos)
    
        Random access operator


    .. cpp:function:: bool end() const
    
        Check if iterator has reached the end


    .. cpp:function:: MeshEntityIterator end_iterator()
    
        Provide a safeguard iterator pointing beyond the end of an iteration
        process, either iterating over the mesh /or incident entities. Added to
        be bit more like STL iterators, since many algorithms rely on a kind of
        beyond iterator.


    .. cpp:function:: void set_end()
    
        Set pos to end position. To create a kind of mesh.end() iterator.


