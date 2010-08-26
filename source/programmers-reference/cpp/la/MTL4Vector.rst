.. Documentation for the header file dolfin/la/MTL4Vector.h

.. _programmers_reference_cpp_la_Mesh:

MTL4Vector.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: implements a minimal backend for MTL4.

.. cpp:class:: MTL4Vector

    *Parent class*
    
        * :cpp:class:`GenericVector`
        
    .. cpp:function:: MTL4Vector()
    
        Create empty vector

    .. cpp:function:: MTL4Vector(const MTL4Vector& x)
    
        Copy constructor

    .. cpp:function:: const mtl4_vector& vec() const
    
        Return const mtl4_vector reference

    .. cpp:function:: explicit MTL4Vector(uint N)
    
        Create vector of size N

    .. cpp:function:: mtl4_vector& vec()
    
        Return mtl4_vector reference

    .. cpp:function:: virtual MTL4Vector* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: virtual const MTL4Vector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const MTL4Vector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: virtual const MTL4Vector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: virtual const MTL4Vector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: virtual const MTL4Vector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: virtual const MTL4Vector& operator= (const MTL4Vector& x)
    
        Assignment operator

    .. cpp:function:: virtual const MTL4Vector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: virtual const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: virtual double inner(const GenericVector& vector) const
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of vector

    .. cpp:function:: virtual double sum() const
    
        Return sum of values of vector

    .. cpp:function:: virtual double* data()
    
        Return pointer to underlying data (non-const version)

    .. cpp:function:: virtual std::pair<uint, uint> local_range() const
    
        Return local ownership range of a vector

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual uint size() const
    
        Return size of vector

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: virtual void add_local(const Array<double>& values)
    
        Add all values to each entry on local process

    .. cpp:function:: virtual void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: virtual void gather(GenericVector& x, const Array<uint>& indices) const
    
        Gather entries into local vector x

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows) const
    
        Get block of values

    .. cpp:function:: virtual void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: virtual void resize(uint N)
    
        Resize vector to size N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: virtual void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual ~MTL4Vector()
    
        Destructor

