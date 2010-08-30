.. Documentation for the header file dolfin/la/PETScVector.h

.. _programmers_reference_cpp_la_petscvector:

PETScVector.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: PETScVector

    *Parent class*
    
        * :cpp:class:`GenericVector,`
        
    This class provides a simple vector class based on PETSc.
    It is a simple wrapper for a PETSc vector pointer (Vec)
    implementing the GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the PETSc Vec pointer using the function vec() and
    use the standard PETSc interface.

    .. cpp:function:: PETScVector(const PETScVector& x)
    
        Copy constructor

    .. cpp:function:: PETScVector(uint N, std::string type="global")
    
        Create vector of size N

    .. cpp:function:: boost::shared_ptr<Vec> vec() const
    
        Return shared_ptr to PETSc Vec object

    .. cpp:function:: const PETScVector& operator= (const PETScVector& x)
    
        Assignment operator

    .. cpp:function:: explicit PETScVector(boost::shared_ptr<Vec> x)
    
        Create vector from given PETSc Vec pointer

    .. cpp:function:: explicit PETScVector(std::string type="global")
    
        Create empty vector

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual PETScVector* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: virtual const PETScVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const PETScVector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: virtual const PETScVector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: virtual const PETScVector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: virtual const PETScVector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: virtual const PETScVector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: virtual double inner(const GenericVector& v) const
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of vector

    .. cpp:function:: virtual double sum() const
    
        Return sum of values of vector

    .. cpp:function:: virtual double sum(const Array<uint>& rows) const
    
        Return sum of selected rows in vector

    .. cpp:function:: virtual std::pair<uint, uint> local_range() const
    
        Return ownership range of a vector

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

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

    .. cpp:function:: virtual void gather(GenericVector& y, const Array<uint>& indices) const
    
        Gather vector entries into a local vector. If local_indices is
        0, then a local index array is created such that the order of
        the values in the return array is the same as the order in
        global_indices.

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows) const
    
        Get block of values (values may live on any process)

    .. cpp:function:: virtual void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: virtual void get_local(double* block, uint m, const uint* rows) const
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: virtual void resize(uint N)
    
        Resize vector ro size N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: virtual void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual ~PETScVector()
    
        Destructor

