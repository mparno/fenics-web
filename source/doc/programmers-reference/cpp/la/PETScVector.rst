.. Documentation for the header file dolfin/la/PETScVector.h

.. _programmers_reference_cpp_la_petscvector:

PETScVector.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: PETScVector

    *Parent class*
    
        * :cpp:class:`GenericVector,`
        
    This class provides a simple vector class based on PETSc.
    It is a simple wrapper for a PETSc vector pointer (Vec)
    implementing the GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the PETSc Vec pointer using the function vec() and
    use the standard PETSc interface.

    .. cpp:function:: explicit PETScVector(std::string type="global")
    
        Create empty vector

    .. cpp:function:: PETScVector(uint N, std::string type="global")
    
        Create vector of size N

    .. cpp:function:: PETScVector(const GenericSparsityPattern& sparsity_pattern)
    
        Create vector

    .. cpp:function:: PETScVector(const PETScVector& x)
    
        Copy constructor

    .. cpp:function:: explicit PETScVector(boost::shared_ptr<Vec> x)
    
        Create vector from given PETSc Vec pointer

    .. cpp:function:: PETScVector* copy() const
    
        Return copy of tensor

    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: void resize(uint N)
    
        Resize vector to global size N

    .. cpp:function:: void resize(std::pair<uint, uint> range)
    
        Resize vector with given ownership range

    .. cpp:function:: void resize(std::pair<uint, uint> range, const std::vector<uint>& ghost_indices)
    
        Resize vector with given ownership range and with ghost values

    .. cpp:function:: uint size() const
    
        Return size of vector

    .. cpp:function:: uint local_size() const
    
        Return local size of vector

    .. cpp:function:: std::pair<uint, uint> local_range() const
    
        Return ownership range of a vector

    .. cpp:function:: bool owns_index(uint i) const
    
        Determine whether global vector index is owned by this process

    .. cpp:function:: void get_local(double* block, uint m, const uint* rows) const
    
        Get block of values (values must all live on the local process)

    .. cpp:function:: void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: void add_local(const Array<double>& values)
    
        Add values to each entry on local process

    .. cpp:function:: void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: void abs()
    
        Replace all entries in the vector by their absolute values

    .. cpp:function:: double inner(const GenericVector& v) const
    
        Return inner product with given vector

    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of vector

    .. cpp:function:: double min() const
    
        Return minimum value of vector

    .. cpp:function:: double max() const
    
        Return maximum value of vector

    .. cpp:function:: double sum() const
    
        Return sum of values of vector

    .. cpp:function:: double sum(const Array<uint>& rows) const
    
        Return sum of selected rows in vector

    .. cpp:function:: const PETScVector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: const PETScVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: const PETScVector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: const PETScVector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: const PETScVector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: const PETScVector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: boost::shared_ptr<Vec> vec() const
    
        Return shared_ptr to PETSc Vec object

    .. cpp:function:: const PETScVector& operator= (const PETScVector& x)
    
        Assignment operator

    .. cpp:function:: void gather(GenericVector& y, const Array<uint>& indices) const
    
        Gather vector entries into a local vector

    .. cpp:function:: void gather(Array<double>& x, const Array<uint>& indices) const
    
        Gather entries into Array x

