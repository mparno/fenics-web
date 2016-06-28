
.. Documentation for the header file dolfin/la/PETScMatrix.h

.. _programmers_reference_cpp_la_petscmatrix:

PETScMatrix.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScMatrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
        * :cpp:class:`PETScBaseMatrix`
        
    This class provides a simple matrix class based on PETSc.
    It is a wrapper for a PETSc matrix pointer (Mat)
    implementing the GenericMatrix interface.
    
    The interface is intentionally simple. For advanced usage,
    access the PETSc Mat pointer using the function mat() and
    use the standard PETSc interface.


    .. cpp:function:: PETScMatrix()
    
        Create empty matrix (on MPI_COMM_WORLD)


    .. cpp:function:: explicit PETScMatrix(MPI_Comm comm)
    
        Create empty matrix


    .. cpp:function:: explicit PETScMatrix(Mat A)
    
        Create a wrapper around a PETSc Mat pointer. The Mat object
        should have been created, e.g. via PETSc MatrCreate.


    .. cpp:function:: PETScMatrix(const PETScMatrix& A)
    
        Copy constructor


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using tensor layout


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::int64_t, std::int64_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: std::size_t nnz() const
    
        Return number of non-zero entries in matrix (collective)


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor. The following values are recognized
        for the mode parameter:
        
          add    - corresponds to PETSc MatAssemblyBegin+End(MAT_FINAL_ASSEMBLY)
          insert - corresponds to PETSc MatAssemblyBegin+End(MAT_FINAL_ASSEMBLY)
          flush  - corresponds to PETSc MatAssemblyBegin+End(MAT_FLUSH_ASSEMBLY)


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<GenericMatrix> copy() const
    
        Return copy of matrix


    .. cpp:function:: void init_vector(GenericVector& z, std::size_t dim) const
    
        Initialize vector z to be compatible with the matrix-vector product
        y = Ax. In the parallel case, both size and layout are
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
    
        Set given rows (global row indices) to identity matrix


    .. cpp:function:: void ident_local(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows (local row indices) to identity matrix


    .. cpp:function:: void get_diagonal(GenericVector& x) const
    
        Get diagonal of a matrix


    .. cpp:function:: void set_diagonal(const GenericVector& x)
    
        Set diagonal of a matrix


    .. cpp:function:: const PETScMatrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const PETScMatrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator


    .. cpp:function:: bool is_symmetric(double tol) const
    
        Test if matrix is symmetric


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: void set_options_prefix(std::string options_prefix)
    
        Sets the prefix used by PETSc when searching the options
        database


    .. cpp:function:: std::string get_options_prefix() const
    
        Returns the prefix used by PETSc when searching the options
        database


    .. cpp:function:: void set_from_options()
    
        Call PETSc function MatSetFromOptions on the PETSc Mat object


    .. cpp:function:: const PETScMatrix& operator= (const PETScMatrix& A)
    
        Assignment operator


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& nullspace)
    
        Attach nullspace to matrix (typically used by Krylov solvers
        when solving singular systems)


    .. cpp:function:: void set_near_nullspace(const VectorSpaceBasis& nullspace)
    
        Attach near nullspace to matrix (used by preconditioners, such
        as smoothed aggregation algerbraic multigrid)


    .. cpp:function:: void binary_dump(std::string file_name) const
    
        Dump matrix to PETSc binary format


