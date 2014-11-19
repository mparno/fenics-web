
.. Documentation for the header file dolfin/la/Matrix.h

.. _programmers_reference_cpp_la_matrix:

Matrix.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Matrix

    *Parent class(es)*
    
        * :cpp:class:`GenericMatrix`
        
    This class provides the default DOLFIN matrix class,
    based on the default DOLFIN linear algebra backend.


    .. cpp:function:: Matrix()
    
        Create empty matrix


    .. cpp:function:: Matrix(const Matrix& A)
    
        Copy constructor


    .. cpp:function:: Matrix(const GenericMatrix& A)
    
        Create a Vector from a GenericVetor


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using tensor layout


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<GenericMatrix> copy() const
    
        Return copy of matrix


    .. cpp:function:: void resize(GenericVector& y, std::size_t dim) const
    
        Resize vector y such that is it compatible with matrix for
        multuplication Ax = b (dim = 0 -> b, dim = 1 -> x) In parallel
        case, size and layout are important.


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


    .. cpp:function:: const Matrix& operator*= (double a)
    
        Multiply matrix by given number


    .. cpp:function:: const Matrix& operator/= (double a)
    
        Divide matrix by given number


    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator


    .. cpp:function:: bool is_symmetric(double tol) const
    
        Test if matrix is symmetric


    .. cpp:function:: boost::tuples::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compressed storage data.
        See GenericMatrix for documentation.


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: const GenericMatrix* instance() const
    
        Return concrete instance / unwrap (const version)


    .. cpp:function:: GenericMatrix* instance()
    
        Return concrete instance / unwrap (non-const version)


    .. cpp:function:: const Matrix& operator= (const Matrix& A)
    
        Assignment operator


