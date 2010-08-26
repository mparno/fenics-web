.. Documentation for the header file dolfin/la/BlockVector.h

.. _programmers_reference_cpp_la_Mesh:

BlockVector.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: BlockVector

    .. cpp:function:: BlockVector(uint n_=0, bool owner=false)
    
        Constructor

    .. cpp:function:: const BlockVector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: const BlockVector& operator+= (const BlockVector& x)
    
        Add given vector

    .. cpp:function:: const BlockVector& operator-= (const BlockVector& x)
    
        Subtract given vector

    .. cpp:function:: const BlockVector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: const BlockVector& operator= (const BlockVector& x)
    
        Assignment operator

    .. cpp:function:: const BlockVector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: const Vector& get(uint i) const
    
        Get functions (const and non-const)

    .. cpp:function:: double inner(const BlockVector& x) const
    
        Return inner product with given vector

    .. cpp:function:: double max() const
    
        Return maximum value of vector

    .. cpp:function:: double min() const
    
        Return minimum value of vector

    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of vector

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: uint size() const
    
        Number of vectors

    .. cpp:function:: virtual BlockVector* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual ~BlockVector()
    
        Destructor

    .. cpp:function:: void axpy(double a, const BlockVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: void set(uint i, Vector& v)
    
        Set function

.. cpp:class:: SubVector

