.. Documentation for the header file dolfin/la/GenericVector.h

.. _programmers_reference_cpp_la_genericvector:

GenericVector.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GenericVector

    *Parent class*
    
        * :cpp:class:`GenericTensor`
        
    This class defines a common interface for vectors.

    .. cpp:function:: void resize(uint rank, const uint* dims)
    
        Resize tensor with given dimensions

    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: GenericVector* copy() const = 0
    
        Return copy of tensor

    .. cpp:function:: uint rank() const
    
        Return tensor rank (number of dimensions)

    .. cpp:function:: uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: std::pair<uint, uint> local_range(uint dim) const
    
        Return local ownership range

    .. cpp:function:: void get(double* block, const uint* num_rows, const uint * const * rows) const
    
        Get block of values

    .. cpp:function:: void set(const double* block, const uint* num_rows, const uint * const * rows)
    
        Set block of values

    .. cpp:function:: void add(const double* block, const uint* num_rows, const uint * const * rows)
    
        Add block of values

    .. cpp:function:: void add(const double* block, const std::vector<const std::vector<uint>* >& rows)
    
        Add block of values

    .. cpp:function:: void add(const double* block, const std::vector<std::vector<uint> >& rows)
    
        Add block of values

    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: void resize(uint N) = 0
    
        Resize vector to global size N

    .. cpp:function:: void resize(std::pair<uint, uint> range) = 0
    
        Resize vector with given ownership range

    .. cpp:function:: void resize(std::pair<uint, uint> range, const std::vector<uint>& ghost_indices) = 0
    
        Resize vector with given ownership range and with ghost values

    .. cpp:function:: uint size() const = 0
    
        Return global size of vector

    .. cpp:function:: uint local_size() const = 0
    
        Return local size of vector

    .. cpp:function:: std::pair<uint, uint> local_range() const = 0
    
        Return local ownership range of a vector

    .. cpp:function:: bool owns_index(uint i) const = 0
    
        Determine whether global vector index is owned by this process

    .. cpp:function:: void get(double* block, uint m, const uint* rows) const
    
        Get block of values (values may live on any process)

    .. cpp:function:: void get_local(double* block, uint m, const uint* rows) const = 0
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: void set(const double* block, uint m, const uint* rows) = 0
    
        Set block of values

    .. cpp:function:: void add(const double* block, uint m, const uint* rows) = 0
    
        Add block of values

    .. cpp:function:: void get_local(Array<double>& values) const = 0
    
        Get all values on local process

    .. cpp:function:: void set_local(const Array<double>& values) = 0
    
        Set all values on local process

    .. cpp:function:: void add_local(const Array<double>& values) = 0
    
        Add values to each entry on local process

    .. cpp:function:: void gather(GenericVector& x, const Array<uint>& indices) const = 0
    
        Gather entries into local vector x

    .. cpp:function:: void gather(Array<double>& x, const Array<uint>& indices) const = 0
    
        Gather entries into Array x

    .. cpp:function:: void axpy(double a, const GenericVector& x) = 0
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: void abs() = 0
    
        Replace all entries in the vector by their absolute values

    .. cpp:function:: double inner(const GenericVector& x) const = 0
    
        Return inner product with given vector

    .. cpp:function:: double norm(std::string norm_type) const = 0
    
        Return norm of vector

    .. cpp:function:: double min() const = 0
    
        Return minimum value of vector

    .. cpp:function:: double max() const = 0
    
        Return maximum value of vector

    .. cpp:function:: double sum() const = 0
    
        Return sum of vector

    .. cpp:function:: double sum(const Array<uint>& rows) const
    
        Return sum of selected rows in vector. Repeated entries only summed once.

    .. cpp:function:: const GenericVector& operator*= (double a) = 0
    
        Multiply vector by given number

    .. cpp:function:: const GenericVector& operator*= (const GenericVector& x) = 0
    
        Multiply vector by another vector pointwise

    .. cpp:function:: const GenericVector& operator/= (double a) = 0
    
        Divide vector by given number

    .. cpp:function:: const GenericVector& operator+= (const GenericVector& x) = 0
    
        Add given vector

    .. cpp:function:: const GenericVector& operator-= (const GenericVector& x) = 0
    
        Subtract given vector

    .. cpp:function:: const GenericVector& operator= (const GenericVector& x) = 0
    
        Assignment operator

    .. cpp:function:: const GenericVector& operator= (double a) = 0
    
        Assignment operator

    .. cpp:function:: const double* data() const
    
        Return pointer to underlying data (const version)

    .. cpp:function:: double* data()
    
        Return pointer to underlying data

    .. cpp:function:: void update_ghost_values()
    
        Update ghost values

    .. cpp:function:: double operator[] (uint i) const
    
        Get value of given entry

    .. cpp:function:: double getitem(uint i) const
    
        Get value of given entry

    .. cpp:function:: void setitem(uint i, double value)
    
        Set given entry to value. apply("insert") should be called before using
        using the object.

