
.. Documentation for the header file dolfin/mesh/MeshEntityIteratorBase.h

.. _programmers_reference_cpp_mesh_meshentityiteratorbase:

MeshEntityIteratorBase.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEntityIteratorBase

    .. cpp:function:: explicit MeshEntityIteratorBase(const Mesh& mesh)
    
        Create iterator for mesh entities over given topological dimension


    .. cpp:function:: MeshEntityIteratorBase(const Mesh& mesh, std::string opt)
    
        Iterator over MeshEntity of dimension dim on mesh, with string option
        to iterate over "regular", "ghost" or "all" entities


    .. cpp:function:: explicit MeshEntityIteratorBase(const MeshEntity& entity)
    
        Create iterator for entities of given dimension connected to given entity


    .. cpp:function:: MeshEntityIteratorBase(const MeshEntityIteratorBase& it)
    
        Copy constructor


    .. cpp:function:: MeshEntityIteratorBase& operator++()
    
        Step to next mesh entity (prefix increment)


    .. cpp:function:: MeshEntityIteratorBase& operator--()
    
        Step to the previous mesh entity (prefix decrease)


    .. cpp:function:: std::size_t pos() const
    
        Return current position


    .. cpp:function:: bool operator==(const MeshEntityIteratorBase & it) const
    
        Comparison operator.


    .. cpp:function:: bool operator!=(const MeshEntityIteratorBase & it) const
    
        Comparison operator


    .. cpp:function:: T& operator*()
    
        Dereference operator


    .. cpp:function:: T* operator->()
    
        Member access operator


    .. cpp:function:: T& operator[] (std::size_t pos)
    
        Random access operator


    .. cpp:function:: bool end() const
    
        Check if iterator has reached the end


    .. cpp:function:: MeshEntityIteratorBase end_iterator()
    
        Provide a safeguard iterator pointing beyond the end of an iteration
        process, either iterating over the mesh /or incident entities. Added to
        be bit more like STL iterators, since many algorithms rely on a kind of
        beyond iterator.


    .. cpp:function:: void set_end()
    
        Set pos to end position. To create a kind of mesh.end() iterator.


