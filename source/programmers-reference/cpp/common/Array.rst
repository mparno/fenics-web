.. Documentation for the header file dolfin/common/Array.h

.. _programmers_reference_cpp_common_array:

Array.h
=======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: T>

    This class provides a simple wrapper for a pointer to an array. A purpose
    of this class is to enable the simple and safe exchange of data between
    C++ and Python.

    .. cpp:function:: Array(const Array& other) : _size(0), x(0)
    
        Copy constructor (arg name need to have a different name that 'x')

    .. cpp:function:: Array(uint N, T* x) : _size(N), x(boost::shared_array<T>(x, NoDeleter<T>()))
    
        Construct array from a pointer. Array will not take ownership.

    .. cpp:function:: Array(uint N, boost::shared_array<T> x) : _size(N), x(x)
    
        Construct array from a shared pointer

    .. cpp:function:: T max() const
    
        Return maximum value of array

    .. cpp:function:: T min() const
    
        Return minimum value of array

    .. cpp:function:: T& operator[] (uint i)
    
        Access value of given entry (non-const version)

    .. cpp:function:: boost::shared_array<T> data()
    
        Return pointer to data (non-const version)

    .. cpp:function:: const Array& operator= (const Array& x)
    
        Assignment operator

    .. cpp:function:: const T& operator[] (uint i) const
    
        Access value of given entry (const version)

    .. cpp:function:: const boost::shared_array<T> data() const
    
        Return pointer to data (const version)

    .. cpp:function:: explicit Array(uint N) : _size(N), x(new T[N])
    
        Create array of size N

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print).
        Note that the Array class is not a subclass of Variable (for
        efficiency) which means that one needs to call str() directly
        instead of using the info() function on Array objects.

    .. cpp:function:: uint size() const
    
        Return size of array

    .. cpp:function:: void resize(uint N)
    
        Resize array to size N. If size changes, contents will be destroyed.

    .. cpp:function:: void update(uint N, T* _x)
    
        Construct array from a pointer. Array will not take ownership.

    .. cpp:function:: void zero()
    
        Zero array

    .. cpp:function:: void zero_eps(double eps=DOLFIN_EPS)
    
        Set entries which meet (abs(x[i]) < eps) to zero

    .. cpp:function:: ~Array()
    
        Destructor

