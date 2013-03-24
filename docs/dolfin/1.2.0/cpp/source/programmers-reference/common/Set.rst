
.. Documentation for the header file dolfin/common/Set.h

.. _programmers_reference_cpp_common_set:

Set.h
=====

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Set

    .. cpp:function:: Set()
    
        Create empty set


    .. cpp:function:: Set(std::vector<T>& x)
    
        Wrap std::vectpr as a set. Contents will be erased.


    .. cpp:function:: Set(const dolfin::Set<T>& x)
    
        Copy constructor


    .. cpp:function:: iterator find(const T& x)
    
        Find entry in set and return an iterator to the entry


    .. cpp:function:: const_iterator find(const T& x) const
    
        Find entry in set and return an iterator to the entry (const)


    .. cpp:function:: bool insert(const T& x)
    
        Insert entry


    .. cpp:function:: void insert(const InputIt first, const InputIt last)
    
        Insert entries


    .. cpp:function:: std::size_t size() const
    
        Set size


    .. cpp:function:: void erase(const T& x)
    
        Erase an entry


    .. cpp:function:: void sort()
    
        Sort set


    .. cpp:function:: void clear()
    
        Clear set


    .. cpp:function:: T operator[](std::size_t n) const
    
        Index the nth entry in the set


    .. cpp:function:: const std::vector<T>& set() const
    
        Return the vector that stores the data in the Set


    .. cpp:function:: std::vector<T>& set()
    
        Return the vector that stores the data in the Set


