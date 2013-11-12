
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


    .. cpp:function:: PETScMatrix(bool use_gpu=false)
    
        Create empty matrix


    .. cpp:function:: explicit PETScMatrix(boost::shared_ptr<Mat> A, bool use_gpu=false)
    
        Create matrix from given PETSc Mat pointer


    .. cpp:function:: PETScMatrix(const PETScMatrix& A)
    
        Copy constructor


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using tensor layout


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor. The following values are recognized
        for the mode parameter:
        
          add    - corresponds to PETSc MatAssemblyBegin+End(MAT_FINAL_ASSEMBLY)
          insert - corresponds to PETSc MatAssemblyBegin+End(MAT_FINAL_ASSEMBLY)
          flush  - corresponds to PETSc MatAssemblyBegin+End(MAT_FLUSH_ASSEMBLY)


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<GenericMatrix> copy() const
    
        Return copy of matrix


    .. cpp:function:: void resize(GenericVector& z, std::size_t dim) const
    
        Resize vector z to be compatible with the matrix-vector product
        y = Ax. In the parallel case, both size and layout are
        important.
        
        *Arguments*
            dim (std::size_t)
                The dimension (axis): dim = 0 --> z = y, dim = 1 --> z = x


    .. cpp:function:: void get(double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols) const
    
        Get block of values


    .. cpp:function:: void set(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Set block of values


    .. cpp:function:: void add(const double* block, std::size_t m, const dolfin::la_index* rows, std::size_t n, const dolfin::la_index* cols)
    
        Add block of values


    .. cpp:function:: void axpy(double a, const GenericMatrix& A, bool same_nonzero_pattern)
    
        Add multiple of given matrix (AXPY operation)


    .. cpp:function:: void getrow(std::size_t row, std::vector<std::size_t>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row


    .. cpp:function:: void setrow(std::size_t row, const std::vector<std::size_t>& columns, const std::vector<double>& values)
    
        Set values for given row


    .. cpp:function:: void zero(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to zero


    .. cpp:function:: void ident(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


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


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix


    .. cpp:function:: const PETScMatrix& operator= (const PETScMatrix& A)
    
        Assignment operator


    .. cpp:function:: void binary_dump(std::string file_name) const
    
        Dump matrix to PETSc binary format


