.. Documentation for the header file dolfin/common/Set.h

.. _programmers_reference_cpp_common_set:

Set.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Set

    This is a std::set like data structure. It is not ordered and it is based
    a std::vector. It can be faster than a std::set for some cases.

    .. cpp:function:: Set(const dolfin::Set<T>& x)
    
        Copy constructor

    .. cpp:function:: Set(std::vector<T>& x)
    
        Wrap std::vectpr as a set. Contents will be erased.

    .. cpp:function:: T operator[](uint n) const
    
        Index the nth entry in the set

    .. cpp:function:: bool insert(const T& x)
    
        Insert entry

    .. cpp:function:: const std::vector<T>& set() const
    
        Return the vector that stores the data in the Set

    .. cpp:function:: const_iterator find(const T& x) const
    
        Find entry in set and return an iterator to the entry (const)

    .. cpp:function:: dolfin::uint size() const
    
        Set size

    .. cpp:function:: iterator find(const T& x)
    
        Find entry in set and return an iterator to the entry

    .. cpp:function:: std::vector<T>& set()
    
        Return the vector that stores the data in the Set

    .. cpp:function:: template<class T> Set()
    
        Create empty set

    .. cpp:function:: void clear()
    
        Clear set

    .. cpp:function:: void erase(const T& x)
    
        Erase an entry

    .. cpp:function:: void resize(uint n)
    
        Resize set

    .. cpp:function:: void sort()
    
        Sort set

