
.. Documentation for the header file dolfin/mesh/SubsetIterator.h

.. _programmers_reference_cpp_mesh_subsetiterator:

SubsetIterator.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SubsetIterator

    A :cpp:class:`SubsetIterator` is similar to a :cpp:class:`MeshEntityIterator` but
    iterates over a specified subset of the range of entities as
    specified by a :cpp:class:`MeshFunction` that labels the entites.


    .. cpp:function:: SubsetIterator(const MeshFunction<std::size_t>& labels, std::size_t label)
    
        Create iterator for given mesh function. The iterator visits
        all entities that match the given label.


    .. cpp:function:: SubsetIterator(const SubsetIterator& subset_iter)
    
        Copy Constructor


    .. cpp:function:: SubsetIterator& operator++()
    
        Step to next mesh entity (prefix increment)


    .. cpp:function:: bool operator==(const SubsetIterator& sub_iter) const
    
        Comparison operator


    .. cpp:function:: bool operator!=(const SubsetIterator & sub_iter) const
    
        Comparison operator


    .. cpp:function:: MeshEntity& operator*()
    
        Dereference operator


    .. cpp:function:: MeshEntity* operator->()
    
        Member access operator


    .. cpp:function:: bool end() const
    
        Check if iterator has reached the end


    .. cpp:function:: void set_end()
    
        Set pos to end position. To create a kind of mesh.end() iterator.


