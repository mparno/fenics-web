
.. Documentation for the header file dolfin/la/EigenMatrix.h

.. _programmers_reference_cpp_la_eigenmatrix:

EigenMatrix.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EigenMatrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
    This class provides a sparse matrix class based on Eigen.  It is
    a simple wrapper for Eigen::SparseMatrix implementing the
    GenericMatrix interface.
    
    The interface is intentionally simple. For advanced usage,
    access the underlying Eigen matrix and use the standard Eigen
    interface which is documented at http://eigen.tuxfamily.org


    .. cpp:function:: EigenMatrix()
    
        Create empty matrix


    .. cpp:function:: EigenMatrix(std::size_t M, std::size_t N)
    
        Create M x N matrix


    .. cpp:function:: EigenMatrix(const EigenMatrix& A)
    
        Copy constructor


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using tenor layout


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::int64_t, std::int64_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: std::size_t nnz() const
    
        Return number of non-zero entries in matrix


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<GenericMatrix> copy() const
    
        Return copy of matrix


    .. cpp:function:: void resize(std::size_t M, std::size_t N)
    
        Resize matrix to M x N


    .. cpp:function:: void init_vector(GenericVector& z, std::size_t dim) const
    
        Initialise vector z to be compatible with the matrix-vector product
        y = Ax.
        
        *Arguments*
            dim (std::size_t)
                The dimension (axis): dim = 0 --> z = y, dim = 1 --> z = x


    .. cpp:function:: void get(double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols) const
    
        Get block of values


    .. cpp:function:: void set(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Set block of values using global indices


    .. cpp:function:: void set_local(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Set block of values using local indices


    .. cpp:function:: void add(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Add block of values using global indices


    .. cpp:function:: void add_local(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Add block of values using local indices


    .. cpp:function:: void axpy(double a, const GenericMatrix& A, bool same_nonzero_pattern)
    
        Add multiple of given matrix (AXPY operation)


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix


    .. cpp:function:: void getrow(std::size_t row, std::vector<std::size_t>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row


    .. cpp:function:: void setrow(std::size_t row_idx, const std::vector<std::size_t>& columns, const std::vector<double>& values)
    
        Set values for given row


    .. cpp:function:: void zero(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (global row indices) to zero


    .. cpp:function:: void zero_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (local row indices) to zero


    .. cpp:function:: void ident(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


    .. cpp:function:: void ident_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const
    
        Matrix-vector product, y = Ax


    .. cpp:function:: void transpmult(const GenericVector& x, GenericVector& y) const
    
        Matrix-vector product, y = A^T x


    .. cpp:function:: void get_diagonal(GenericVector& x) const
    
        Get diagonal of a matrix


    .. cpp:function:: void set_diagonal(const GenericVector& x)
    
        Set diagonal of a matrix


    .. cpp:function:: const EigenMatrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const EigenMatrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator


    .. cpp:function:: std::tuple<const int*, const int*, const double*, std::size_t> data() const
    
        Return pointers to underlying compressed storage data See
        GenericMatrix for documentation.


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: const eigen_matrix_type& mat() const
    
        Return reference to Eigen matrix (const version)


    .. cpp:function:: eigen_matrix_type& mat()
    
        Return reference to Eigen matrix (non-const version)


    .. cpp:function:: void compress()
    
        Compress matrix (eliminate all zeros from a sparse matrix)


    .. cpp:function:: double operator() (dolfin::la_index i, dolfin::la_index j) const
    
        Access value of given entry


    .. cpp:function:: const EigenMatrix& operator= (const EigenMatrix& A)
    
        Assignment operator


