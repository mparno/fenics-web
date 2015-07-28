
.. Documentation for the header file dolfin/la/STLMatrix.h

.. _programmers_reference_cpp_la_stlmatrix:

STLMatrix.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: STLMatrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
    Simple STL-based implementation of the GenericMatrix interface.
    The sparse matrix is stored as a pair of std::vector of
    std::vector, one for the columns and one for the values.
    
    Historically, this class has undergone a number of different
    incarnations, based on various combinations of std::vector,
    std::set and std::map. The current implementation has proven to
    be the fastest.


    .. cpp:function:: STLMatrix(std::size_t primary_dim=0)
    
        Create empty matrix


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        --- Implementation of the GenericTensor interface ---
        Initialize zero tensor using sparsity pattern


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local ownership range


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


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix


    .. cpp:function:: void getrow(std::size_t row, std::vector<std::size_t>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row


    .. cpp:function:: void setrow(std::size_t row, const std::vector<std::size_t>& columns, const std::vector<double>& values)
    
        Set values for given row


    .. cpp:function:: void zero(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (global row indices) to zero


    .. cpp:function:: void zero_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (local row indices) to zero


    .. cpp:function:: void ident(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


    .. cpp:function:: void ident_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


    .. cpp:function:: void get_diagonal(GenericVector& x) const
    
        Get diagonal of a matrix


    .. cpp:function:: void set_diagonal(const GenericVector& x)
    
        Set diagonal of a matrix


    .. cpp:function:: const STLMatrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const STLMatrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        --- Specialized matrix functions ---
        Return linear algebra backend factory


    .. cpp:function:: std::size_t block_size() const
    
        --- STLMatrix interface ---
        Return matrix block size


    .. cpp:function:: void clear()
    
        Clear matrix. Destroys data and sparse layout


    .. cpp:function:: void csr(std::vector<double>& vals, std::vector<T>& cols, std::vector<T>& row_ptr, std::vector<T>& local_to_global_row, bool block, bool symmetric) const
    
        Return matrix in CSR format


    .. cpp:function:: void csc(std::vector<double>& vals, std::vector<T>& rows, std::vector<T>& col_ptr, std::vector<T>& local_to_global_col, bool block, bool symmetric) const
    
        Return matrix in CSC format


    .. cpp:function:: std::size_t nnz() const
    
        Return number of global non-zero entries


    .. cpp:function:: std::size_t local_nnz() const
    
        Return number of local non-zero entries


    .. cpp:function:: void compressed_storage(std::vector<double>& vals, std::vector<T>& rows, std::vector<T>& col_ptr, std::vector<T>& local_to_global_col, bool block, bool symmetric) const
    
        Return matrix in compressed format


