
.. Documentation for the header file dolfin/la/PETScVector.h

.. _programmers_reference_cpp_la_petscvector:

PETScVector.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScVector

    *Parent class(es)*
    
        * :cpp:class:`GenericVector`
        
        * :cpp:class:`PETScObject`
        
    This class provides a simple vector class based on PETSc.
    It is a simple wrapper for a PETSc vector pointer (Vec)
    implementing the GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the PETSc Vec pointer using the function vec() and
    use the standard PETSc interface.


    .. cpp:function:: PETScVector()
    
        Create empty vector (on MPI_COMM_WORLD)


    .. cpp:function:: explicit PETScVector(MPI_Comm comm)
    
        Create empty vector on an MPI communicator


    .. cpp:function:: PETScVector(MPI_Comm comm, std::size_t N)
    
        Create vector of size N


    .. cpp:function:: explicit PETScVector(const SparsityPattern& sparsity_pattern)
    
        Create vector


    .. cpp:function:: PETScVector(const PETScVector& x)
    
        Copy constructor


    .. cpp:function:: explicit PETScVector(Vec x)
    
        Create vector wrapper of PETSc Vec pointer. The reference
        counter of the Vec will be increased, and decreased upon
        destruction of this object.


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<GenericVector> copy() const
    
        Return copy of vector


    .. cpp:function:: void init(MPI_Comm comm, std::size_t N)
    
        Initialize vector to global size N


    .. cpp:function:: void init(MPI_Comm comm, std::pair<std::size_t, std::size_t> range)
    
        Initialize vector with given ownership range


    .. cpp:function:: void init(MPI_Comm comm, std::pair<std::size_t, std::size_t> range, const std::vector<std::size_t>& local_to_global_map, const std::vector<la_index>& ghost_indices)
    
        Initialize vector with given ownership range and with ghost
        values


    .. cpp:function:: bool empty() const
    
        Return true if vector is empty


    .. cpp:function:: std::size_t size() const
    
        Return size of vector


    .. cpp:function:: std::size_t local_size() const
    
        Return local size of vector


    .. cpp:function:: std::pair<std::int64_t, std::int64_t> local_range() const
    
        Return ownership range of a vector


    .. cpp:function:: bool owns_index(std::size_t i) const
    
        Determine whether global vector index is owned by this process


    .. cpp:function:: void get(double* block, std::size_t m, const dolfin::la_index* rows) const
    
        Get block of values using global indices (all values must be
        owned by local process, ghosts cannot be accessed)


    .. cpp:function:: void get_local(double* block, std::size_t m, const dolfin::la_index* rows) const
    
        Get block of values using local indices


    .. cpp:function:: void set(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Set block of values using global indices


    .. cpp:function:: void set_local(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Set block of values using local indices


    .. cpp:function:: void add(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Add block of values using global indices


    .. cpp:function:: void add_local(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Add block of values using local indices


    .. cpp:function:: void get_local(std::vector<double>& values) const
    
        Get all values on local process


    .. cpp:function:: void set_local(const std::vector<double>& values)
    
        Set all values on local process


    .. cpp:function:: void add_local(const Array<double>& values)
    
        Add values to each entry on local process


    .. cpp:function:: void gather(GenericVector& y, const std::vector<dolfin::la_index>& indices) const
    
        Gather vector entries into a local vector


    .. cpp:function:: void gather(std::vector<double>& x, const std::vector<dolfin::la_index>& indices) const
    
        Gather entries into x


    .. cpp:function:: void gather_on_zero(std::vector<double>& x) const
    
        Gather all entries into x on process 0


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


    .. cpp:function:: double sum(const Array<std::size_t>& rows) const
    
        Return sum of selected rows in vector


    .. cpp:function:: const PETScVector& operator*= (double a)
    
        Multiply vector by given number


    .. cpp:function:: const PETScVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise


    .. cpp:function:: const PETScVector& operator/= (double a)
    
        Divide vector by given number


    .. cpp:function:: const PETScVector& operator+= (const GenericVector& x)
    
        Add given vector


    .. cpp:function:: const PETScVector& operator+= (double a)
    
        Add number to all components of a vector


    .. cpp:function:: const PETScVector& operator-= (const GenericVector& x)
    
        Subtract given vector


    .. cpp:function:: const PETScVector& operator-= (double a)
    
        Subtract number from all components of a vector


    .. cpp:function:: const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator


    .. cpp:function:: const PETScVector& operator= (double a)
    
        Assignment operator


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: void set_options_prefix(std::string options_prefix)
    
        Sets the prefix used by PETSc when searching the options
        database


    .. cpp:function:: std::string get_options_prefix() const
    
        Returns the prefix used by PETSc when searching the options
        database


    .. cpp:function:: void set_from_options()
    
        Call PETSc function VecSetFromOptions on the underlying Vec
        object


    .. cpp:function:: Vec vec() const
    
        Return pointer to PETSc Vec object


    .. cpp:function:: const PETScVector& operator= (const PETScVector& x)
    
        Assignment operator


    .. cpp:function:: void reset(Vec vec)
    
        Switch underlying PETSc object. Intended for internal library
        usage.


