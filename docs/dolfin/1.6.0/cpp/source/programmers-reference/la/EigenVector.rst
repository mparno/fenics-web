
.. Documentation for the header file dolfin/la/EigenVector.h

.. _programmers_reference_cpp_la_eigenvector:

EigenVector.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EigenVector

    *Parent class(es)*
    
        * :cpp:class:`GenericVector`
        
    This class provides a simple vector class based on Eigen.
    It is a simple wrapper for a Eigen vector implementing the
    GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the underlying Eigen vector and use the standard Eigen
    interface which is documented at http://eigen.tuxfamily.org


    .. cpp:function:: EigenVector()
    
        Create empty vector


    .. cpp:function:: explicit EigenVector(std::size_t N)
    
        Create vector of size N


    .. cpp:function:: EigenVector(const EigenVector& x)
    
        Copy constructor


    .. cpp:function:: explicit EigenVector(std::shared_ptr<Eigen::VectorXd> x)
    
        Construct vector from an Eigen shared_ptr


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<GenericVector> copy() const
    
        Create copy of tensor


    .. cpp:function:: void init(MPI_Comm comm, std::size_t N)
    
        Initialize vector to size N


    .. cpp:function:: void init(MPI_Comm comm, std::pair<std::size_t, std::size_t> range)
    
        Resize vector with given ownership range


    .. cpp:function:: void init(MPI_Comm comm, std::pair<std::size_t, std::size_t> range, const std::vector<std::size_t>& local_to_global_map, const std::vector<la_index>& ghost_indices)
    
        Resize vector with given ownership range and with ghost values


    .. cpp:function:: bool empty() const
    
        Return true if vector is empty


    .. cpp:function:: std::size_t size() const
    
        Return true if vector is empty


    .. cpp:function:: std::size_t local_size() const
    
        Return local size of vector


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range() const
    
        Return local ownership range of a vector


    .. cpp:function:: bool owns_index(std::size_t i) const
    
        Determine whether global vector index is owned by this process


    .. cpp:function:: void get(double* block, std::size_t m, const dolfin::la_index* rows) const
    
        Get block of values using global indices


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


    .. cpp:function:: void gather(GenericVector& x, const std::vector<dolfin::la_index>& indices) const
    
        Gather entries into local vector x


    .. cpp:function:: void gather(std::vector<double>& x, const std::vector<dolfin::la_index>& indices) const
    
        Gather entries into x


    .. cpp:function:: void gather_on_zero(std::vector<double>& x) const
    
        Gather all entries into x on process 0


    .. cpp:function:: void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)


    .. cpp:function:: void abs()
    
        Replace all entries in the vector by their absolute values


    .. cpp:function:: double inner(const GenericVector& x) const
    
        Return inner product with given vector


    .. cpp:function:: double norm(std::string norm_type) const
    
        Compute norm of vector


    .. cpp:function:: double min() const
    
        Return minimum value of vector


    .. cpp:function:: double max() const
    
        Return maximum value of vector


    .. cpp:function:: double sum() const
    
        Return sum of values of vector


    .. cpp:function:: double sum(const Array<std::size_t>& rows) const
    
        Return sum of selected rows in vector. Repeated entries are
        only summed once.


    .. cpp:function:: const EigenVector& operator*= (double a)
    
        Multiply vector by given number


    .. cpp:function:: const EigenVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise


    .. cpp:function:: const EigenVector& operator/= (double a)
    
        Divide vector by given number


    .. cpp:function:: const EigenVector& operator+= (const GenericVector& x)
    
        Add given vector


    .. cpp:function:: const EigenVector& operator+= (double a)
    
        Add number to all components of a vector


    .. cpp:function:: const EigenVector& operator-= (const GenericVector& x)
    
        Subtract given vector


    .. cpp:function:: const EigenVector& operator-= (double a)
    
        Subtract number from all components of a vector


    .. cpp:function:: const GenericVector& operator= (const GenericVector& x)
    
        Assignment operator


    .. cpp:function:: const EigenVector& operator= (double a)
    
        Assignment operator


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: void resize(std::size_t N)
    
        Resize vector to size N


    .. cpp:function:: const Eigen::VectorXd& vec() const
    
        Return reference to Eigen vector (const version)


    .. cpp:function:: Eigen::VectorXd& vec()
    
        Return reference to Eigen vector (non-const version)


    .. cpp:function:: double operator[] (dolfin::la_index i) const
    
        Access value of given entry (const version)


    .. cpp:function:: double& operator[] (dolfin::la_index i)
    
        Access value of given entry (non-const version)


    .. cpp:function:: const EigenVector& operator= (const EigenVector& x)
    
        Assignment operator


    .. cpp:function:: double* data()
    
        Return pointer to underlying data


    .. cpp:function:: const double* data() const
    
        Return pointer to underlying data (const version)


