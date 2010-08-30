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

    .. cpp:function:: Vector()
    
        Create empty vector

    .. cpp:function:: Vector(const GenericVector& x)
    
        Create a Vector from a GenericVetor

    .. cpp:function:: Vector(const Vector& x)
    
        Copy constructor

    .. cpp:function:: const Vector& operator= (const Vector& x)
    
        Assignment operator

    .. cpp:function:: const Vector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: explicit Vector(uint N)
    
        Create vector of size N

    .. cpp:function:: virtual GenericVector* instance()
    
        Return concrete instance / unwrap (non-const version)

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual Vector* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: virtual const GenericVector* instance() const
    
        Return concrete instance / unwrap (const version)

    .. cpp:function:: virtual const Vector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const Vector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: virtual const Vector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: virtual const Vector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: virtual const Vector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: virtual const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: virtual double inner(const GenericVector& x) const
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std
    
        Return norm of vector

    .. cpp:function:: virtual double sum() const
    
        Return sum of values of vector

    .. cpp:function:: virtual double* data()
    
        Return pointer to underlying data

    .. cpp:function:: virtual std
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual std
    
        Return local ownership range of a vector

    .. cpp:function:: virtual uint size() const
    
        Return size of vector

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: virtual void add_local(const Array<double>& values)
    
        Add values to each entry on local process

    .. cpp:function:: virtual void apply(std
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: virtual void gather(GenericVector& x, const Array<uint>& indices) const
    
        Gather entries into local vector x

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows) const
    
        Get block of values

    .. cpp:function:: virtual void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: virtual void get_local(double* block, uint m, const uint* rows) const
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: virtual void resize(uint N)
    
        Resize vector to size N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: virtual void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual ~Vector()
    
        Destructor

