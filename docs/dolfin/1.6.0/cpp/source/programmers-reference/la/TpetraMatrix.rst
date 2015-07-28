
.. Documentation for the header file dolfin/la/TpetraMatrix.h

.. _programmers_reference_cpp_la_tpetramatrix:

TpetraMatrix.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TpetraMatrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
    This class provides a simple matrix class based on Tpetra.  It
    is a wrapper for a Tpetra matrix pointer
    (Teuchos::RCP<matrix_type>) implementing the GenericMatrix
    interface.
    
    The interface is intentionally simple. For advanced usage,
    access the Tpetra::RCP<matrix_type> pointer using the function
    mat() and use the standard Tpetra interface.


    .. cpp:function:: TpetraMatrix()
    
        Create empty matrix


    .. cpp:function:: explicit TpetraMatrix(Teuchos::RCP<matrix_type> A)
    
        Create a wrapper around a Teuchos::RCP<matrix_type> pointer


    .. cpp:function:: TpetraMatrix(const TpetraMatrix& A)
    
        Copy constructor


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using tensor layout


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor. The mode parameter is ignored.


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<GenericMatrix> copy() const
    
        Return copy of matrix


    .. cpp:function:: void init_vector(GenericVector& z, std::size_t dim) const
    
        Initialize vector z to be compatible with the matrix-vector
        product y = Ax. In the parallel case, both size and layout are
        important.
        
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


    .. cpp:function:: void getrow(std::size_t row, std::vector<std::size_t>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row


    .. cpp:function:: void setrow(std::size_t row, const std::vector<std::size_t>& columns, const std::vector<double>& values)
    
        Set values for given row


    .. cpp:function:: void zero(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (global row indices) to zero


    .. cpp:function:: void zero_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (local row indices) to zero


    .. cpp:function:: void ident(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (global row indices) to identity matrix


    .. cpp:function:: void ident_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (local row indices) to identity matrix


    .. cpp:function:: void get_diagonal(GenericVector& x) const
    
        Get diagonal of a matrix


    .. cpp:function:: void set_diagonal(const GenericVector& x)
    
        Set diagonal of a matrix


    .. cpp:function:: const TpetraMatrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const TpetraMatrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator


    .. cpp:function:: bool is_symmetric(double tol) const
    
        Test if matrix is symmetric


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix


    .. cpp:function:: const TpetraMatrix& operator= (const TpetraMatrix& A)
    
        Assignment operator


