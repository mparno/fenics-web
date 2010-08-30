.. Documentation for the header file dolfin/la/Vector.h

.. _programmers_reference_cpp_la_vector:

Vector.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Vector

    *Parent class*
    
        * :cpp:class:`GenericVector`
        
    This class provides the default DOLFIN vector class,
    based on the default DOLFIN linear algebra backend.

    .. cpp:function:: GenericVector* instance()
    
        Return concrete instance / unwrap (non-const version)

    .. cpp:function:: LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: Vector()
    
        Create empty vector

    .. cpp:function:: Vector(const GenericVector& x)
    
        Create a Vector from a GenericVetor

    .. cpp:function:: Vector(const Vector& x)
    
        Copy constructor

    .. cpp:function:: Vector* copy() const
    
        Return copy of tensor

    .. cpp:function:: const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: const GenericVector* instance() const
    
        Return concrete instance / unwrap (const version)

    .. cpp:function:: const Vector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: const Vector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: const Vector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: const Vector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: const Vector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: const Vector& operator= (const Vector& x)
    
        Assignment operator

    .. cpp:function:: const Vector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: double inner(const GenericVector& x) const
    
        Return inner product with given vector

    .. cpp:function:: double max() const
    
        Return maximum value of vector

    .. cpp:function:: double min() const
    
        Return minimum value of vector

    .. cpp:function:: double norm(std
    
        Return norm of vector

    .. cpp:function:: double sum() const
    
        Return sum of values of vector

    .. cpp:function:: double* data()
    
        Return pointer to underlying data

    .. cpp:function:: explicit Vector(uint N)
    
        Create vector of size N

    .. cpp:function:: std
    
        Return informal string representation (pretty-print)

    .. cpp:function:: std
    
        Return local ownership range of a vector

    .. cpp:function:: uint size() const
    
        Return size of vector

    .. cpp:function:: void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: void add_local(const Array<double>& values)
    
        Add values to each entry on local process

    .. cpp:function:: void apply(std
    
        Finalize assembly of tensor

    .. cpp:function:: void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: void gather(GenericVector& x, const Array<uint>& indices) const
    
        Gather entries into local vector x

    .. cpp:function:: void get(double* block, uint m, const uint* rows) const
    
        Get block of values

    .. cpp:function:: void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: void get_local(double* block, uint m, const uint* rows) const
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: void resize(uint N)
    
        Resize vector to size N

    .. cpp:function:: void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: ~Vector()
    
        Destructor

