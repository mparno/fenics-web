
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


    .. cpp:function:: IndexSet(std::size_t size)
    
        Create index set of given size


    .. cpp:function:: bool empty() const
    
        Return true if set is empty


    .. cpp:function:: std::size_t size() const
    
        Return size of set


    .. cpp:function:: bool has_index(std::size_t index) const
    
        Check whether index is in set


    .. cpp:function:: std::size_t find(std::size_t index) const
    
        Return position (if any) for given index


    .. cpp:function:: std::size_t& operator[] (std::size_t i)
    
        Return given index


    .. cpp:function:: const std::size_t& operator[] (std::size_t i) const
    
        Return given index (const version)


    .. cpp:function:: void insert(std::size_t index)
    
        Insert index into set


    .. cpp:function:: void fill()
    
        Fill index set with indices 0, 1, 2, ..., size - 1


    .. cpp:function:: void clear()
    
        Clear set


