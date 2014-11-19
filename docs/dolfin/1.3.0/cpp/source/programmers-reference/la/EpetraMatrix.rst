
.. Documentation for the header file dolfin/la/EpetraMatrix.h

.. _programmers_reference_cpp_la_epetramatrix:

EpetraMatrix.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EpetraMatrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
    This class provides a simple matrix class based on Epetra.
    It is a simple wrapper for an Epetra matrix object (Epetra_FECrsMatrix)
    implementing the GenericMatrix interface.
    
    The interface is intentionally simple. For advanced usage,
    access the Epetra_FECrsMatrix object using the function mat() and
    use the standard Epetra interface.


    .. cpp:function:: EpetraMatrix()
    
        Create empty matrix


    .. cpp:function:: EpetraMatrix(const EpetraMatrix& A)
    
        Copy constuctor


    .. cpp:function:: explicit EpetraMatrix(Teuchos::RCP<Epetra_FECrsMatrix> A)
    
        Create matrix from given Epetra_FECrsMatrix pointer


    .. cpp:function:: explicit EpetraMatrix(boost::shared_ptr<Epetra_FECrsMatrix> A)
    
        Create matrix from given Epetra_FECrsMatrix pointer


    .. cpp:function:: explicit EpetraMatrix(const Epetra_CrsGraph& graph)
    
        Create matrix from given Epetra_CrsGraph


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
        
          add    - corresponding to Epetra GlobalAssemble(Add)
          insert - corresponding to Epetra GlobalAssemble(Insert)


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


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix


    .. cpp:function:: void getrow(std::size_t row, std::vector<std::size_t>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row


    .. cpp:function:: void setrow(std::size_t row, const std::vector<std::size_t>& columns, const std::vector<double>& values)
    
        Set values for given row


    .. cpp:function:: void zero(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to zero


    .. cpp:function:: void ident(std::size_t m, const dolfin::la_index* rows)
    
        Set given rows to identity matrix


    .. cpp:function:: const EpetraMatrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const EpetraMatrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& x)
    
        Assignment operator


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: boost::shared_ptr<Epetra_FECrsMatrix> mat() const
    
        Return Epetra_FECrsMatrix pointer


    .. cpp:function:: const EpetraMatrix& operator= (const EpetraMatrix& x)
    
        Assignment operator


