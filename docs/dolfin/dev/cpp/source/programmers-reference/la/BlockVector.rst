
.. Documentation for the header file dolfin/la/BlockVector.h

.. _programmers_reference_cpp_la_blockvector:

BlockVector.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BlockVector

    .. cpp:function:: BlockVector(std::size_t n = 0)
    
        Constructor


    .. cpp:function:: BlockVector* copy() const
    
        Return copy of tensor


    .. cpp:function:: void set_block(std::size_t i, boost::shared_ptr<GenericVector> v)
    
        Set function


    .. cpp:function:: const boost::shared_ptr<GenericVector> get_block(std::size_t i) const
    
        Get sub-vector (const)


    .. cpp:function:: boost::shared_ptr<GenericVector> get_block(std::size_t)
    
        Get sub-vector (non-const)


    .. cpp:function:: void axpy(double a, const BlockVector& x)
    
        Add multiple of given vector (AXPY operation)


    .. cpp:function:: double inner(const BlockVector& x) const
    
        Return inner product with given vector


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of vector


    .. cpp:function:: double min() const
    
        Return minimum value of vector


    .. cpp:function:: double max() const
    
        Return maximum value of vector


    .. cpp:function:: const BlockVector& operator*= (double a)
    
        Multiply vector by given number


    .. cpp:function:: const BlockVector& operator/= (double a)
    
        Divide vector by given number


    .. cpp:function:: const BlockVector& operator+= (const BlockVector& x)
    
        Add given vector


    .. cpp:function:: const BlockVector& operator-= (const BlockVector& x)
    
        Subtract given vector


    .. cpp:function:: const BlockVector& operator= (const BlockVector& x)
    
        Assignment operator


    .. cpp:function:: const BlockVector& operator= (double a)
    
        Assignment operator


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t size() const
    
        Number of vectors


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


