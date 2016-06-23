
.. Documentation for the header file dolfin/common/ArrayView.h

.. _programmers_reference_cpp_common_arrayview:

ArrayView.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ArrayView

    This class provides a wrapper for a pointer to an array. It
    never owns the data, and will not be valid if the underlying
    data goes out-of-scope.


    .. cpp:function:: ArrayView()
    
        Constructor


    .. cpp:function:: ArrayView(std::size_t N, T* x)
    
        Construct array from a pointer. Array does not take ownership.


    .. cpp:function:: explicit ArrayView(V& v)
    
        Construct array from a container with the the data() and
        size() functions


    .. cpp:function:: ArrayView(const ArrayView& x)
    
        Copy constructor


    .. cpp:function:: void set(std::size_t N, T* x)
    
        Update object to point to new data


    .. cpp:function:: void set(V& v)
    
        Update object to point to new container


    .. cpp:function:: std::size_t size() const
    
        Return size of array


    .. cpp:function:: bool empty() const
    
        Test if array view is empty


    .. cpp:function:: const T& operator[] (std::size_t i) const
    
        Access value of given entry (const version)


    .. cpp:function:: T& operator[] (std::size_t i)
    
        Access value of given entry (non-const version)


    .. cpp:function:: const T* data() const
    
        Return pointer to data (const version)


    .. cpp:function:: T* data()
    
        Return pointer to data (non-const version)


