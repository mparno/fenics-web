.. Documentation for the header file dolfin/la/uBLASVector.h

.. _programmers_reference_cpp_la_Mesh:

uBLASVector.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: uBLASVector

    *Parent class*
    
        * :cpp:class:`GenericVector`
        
        This class provides a simple vector class based on uBLAS.
        It is a simple wrapper for a uBLAS vector implementing the
        GenericVector interface.
        
        The interface is intentionally simple. For advanced usage,
        access the underlying uBLAS vector and use the standard
        uBLAS interface which is documented at
        http://www.boost.org/libs/numeric/ublas/doc/index.htm.

    .. cpp:function:: const uBLASVector& operator= (const uBLASVector& x)
    
        Assignment operator

    .. cpp:function:: const ublas_vector& vec() const
    
        Return reference to uBLAS vector (const version)

    .. cpp:function:: double& operator[] (uint i)
    
        Access value of given entry (non-const version)

    .. cpp:function:: explicit uBLASVector(boost::shared_ptr<ublas_vector> x)
    
        Construct vector from a ublas_vector

    .. cpp:function:: explicit uBLASVector(uint N)
    
        Create vector of size N

    .. cpp:function:: uBLASVector()
    
        Create empty vector

    .. cpp:function:: uBLASVector(const uBLASVector& x)
    
        Copy constructor

    .. cpp:function:: ublas_vector& vec()
    
        Return reference to uBLAS vector (non-const version)

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: virtual const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: virtual const uBLASVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const uBLASVector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: virtual const uBLASVector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: virtual const uBLASVector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: virtual const uBLASVector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: virtual const uBLASVector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: virtual double inner(const GenericVector& x) const
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Compute norm of vector

    .. cpp:function:: virtual double operator[] (uint i) const
    
        Access value of given entry (const version)

    .. cpp:function:: virtual double sum() const
    
        Return sum of values of vector

    .. cpp:function:: virtual double* data()
    
        Return pointer to underlying data

    .. cpp:function:: virtual std::pair<uint, uint> local_range() const
    
        Return local ownership range of a vector

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual uBLASVector* copy() const
    
        Create copy of tensor

    .. cpp:function:: virtual uint size() const
    
        Return size of vector

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: virtual void add_local(const Array<double>& values)
    
        Add values to each entry on local process

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

    .. cpp:function:: virtual ~uBLASVector()
    
        Destructor

