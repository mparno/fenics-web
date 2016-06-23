
.. Documentation for the header file dolfin/common/Array.h

.. _programmers_reference_cpp_common_array:

Array.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Array

    This class provides a simple wrapper for a pointer to an array. A
    purpose of this class is to enable the simple and safe exchange
    of data between C++ and Python.


    .. cpp:function:: explicit Array(std::size_t N)
    
        Create array of size N. Array has ownership.


    .. cpp:function:: Array(std::size_t N, T* x)
    
        Construct array from a pointer. Array does not take ownership.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print).
        Note that the Array class is not a subclass of Variable (for
        efficiency) which means that one needs to call str() directly
        instead of using the info() function on Array objects.


    .. cpp:function:: std::size_t size() const
    
        Return size of array


    .. cpp:function:: const T& operator[] (std::size_t i) const
    
        Access value of given entry (const version)


    .. cpp:function:: T& operator[] (std::size_t i)
    
        Access value of given entry (non-const version)


    .. cpp:function:: const T* data() const
    
        Return pointer to data (const version)


    .. cpp:function:: T* data()
    
        Return pointer to data (non-const version)


    .. cpp:function:: Array(const Array& other) /* leave body undefined */
    
        Disable copy construction, to avoid unanticipated sharing or
        copying of data. This means that an Array must always be passed as
        reference, or as a (possibly shared) pointer.


