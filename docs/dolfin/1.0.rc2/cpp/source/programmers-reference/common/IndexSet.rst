
.. Documentation for the header file dolfin/common/IndexSet.h

.. _programmers_reference_cpp_common_indexset:

IndexSet.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: IndexSet

    This class provides an efficient data structure for index sets.
    The cost of checking whether a given index is in the set is O(1)
    and very very fast (optimal) at the cost of extra storage.


    .. cpp:function:: IndexSet(uint size)
    
        Create index set of given size


    .. cpp:function:: uint size() const
    
        Return size of set


    .. cpp:function:: bool has_index(uint index) const
    
        Check whether index is in set


    .. cpp:function:: uint find(uint index) const
    
        Return position (if any) for given index


    .. cpp:function:: uint& operator[] (uint i)
    
        Return given index


    .. cpp:function:: const uint& operator[] (uint i) const
    
        Return given index (const version)


    .. cpp:function:: void insert(uint index)
    
        Insert index into set


    .. cpp:function:: void fill()
    
        Fill index set with indices 0, 1, 2, ..., size - 1


    .. cpp:function:: void clear()
    
        Clear set


