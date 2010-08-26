.. Documentation for the header file dolfin/la/GenericVector.h

.. _programmers_reference_cpp_la_genericvector:

GenericVector.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: GenericVector

    *Parent class*
    
        * :cpp:class:`GenericTensor`
        
        This class defines a common interface for vectors.

    .. cpp:function:: inline void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: virtual GenericVector* copy() const = 0
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericVector& operator*= (const GenericVector& x) = 0
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const GenericVector& operator*= (double a) = 0
    
        Multiply vector by given number

    .. cpp:function:: virtual const GenericVector& operator+= (const GenericVector& x) = 0
    
        Add given vector

    .. cpp:function:: virtual const GenericVector& operator-= (const GenericVector& x) = 0
    
        Subtract given vector

    .. cpp:function:: virtual const GenericVector& operator/= (double a) = 0
    
        Divide vector by given number

    .. cpp:function:: virtual const GenericVector& operator= (const GenericVector& x) = 0
    
        Assignment operator

    .. cpp:function:: virtual const GenericVector& operator= (double a) = 0
    
        Assignment operator

    .. cpp:function:: virtual const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: virtual double getitem(uint i) const
    
        Get value of given entry

    .. cpp:function:: virtual double inner(const GenericVector& x) const = 0
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const = 0
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const = 0
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std::string norm_type) const = 0
    
        Return norm of vector

    .. cpp:function:: virtual double operator[] (uint i) const
    
        Get value of given entry

    .. cpp:function:: virtual double sum() const = 0
    
        Return sum of vector

    .. cpp:function:: virtual double sum(const Array<uint>& rows) const
    
        Return sum of selected rows in vector. Repeated entries only summed once.

    .. cpp:function:: virtual double* data()
    
        Return pointer to underlying data

    .. cpp:function:: virtual std::pair<uint, uint> local_range() const = 0
    
        Return local ownership range of a vector

    .. cpp:function:: virtual std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual uint local_size() const
    
        Return local size of vector

    .. cpp:function:: virtual uint rank() const
    
        Return tensor rank (number of dimensions)

    .. cpp:function:: virtual uint size() const = 0
    
        Return global size of vector

    .. cpp:function:: virtual uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: virtual void add(const double* block, const uint* num_rows,
                                       const uint * const * rows)
    
        Add block of values

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows) = 0
    
        Add block of values

    .. cpp:function:: virtual void add_local(const Array<double>& values) = 0
    
        Add values to each entry on local process

    .. cpp:function:: virtual void apply(std::string mode) = 0
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericVector& x) = 0
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: virtual void gather(GenericVector& x, const Array<uint>& indices) const = 0
    
        Gather entries into local vector x

    .. cpp:function:: virtual void get(double* block, const uint* num_rows,
                                       const uint * const * rows) const
    
        Get block of values

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows) const = 0
    
        Get block of values (values may live on any process)

    .. cpp:function:: virtual void get_local(Array<double>& values) const = 0
    
        Get all values on local process

    .. cpp:function:: virtual void get_local(double* block, uint m, const uint* rows) const
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: virtual void resize(uint N) = 0
    
        Resize vector to size N

    .. cpp:function:: virtual void resize(uint rank, const uint* dims)
    
        Resize tensor with given dimensions

    .. cpp:function:: virtual void set(const double* block, const uint* num_rows,
                                       const uint * const * rows)
    
        Set block of values

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows) = 0
    
        Set block of values

    .. cpp:function:: virtual void set_local(const Array<double>& values) = 0
    
        Set all values on local process

    .. cpp:function:: virtual void setitem(uint i, double value)
    
        Set given entry to value. apply("insert") should be called before using
        using the object.

    .. cpp:function:: virtual void zero() = 0
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual ~GenericVector()
    
        Destructor

