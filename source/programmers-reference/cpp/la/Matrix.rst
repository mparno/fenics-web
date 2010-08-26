.. Documentation for the header file dolfin/la/Matrix.h

.. _programmers_reference_cpp_la_matrix:

Matrix.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Matrix

    *Parent class*
    
        * :cpp:class:`GenericMatrix`
        
        This class provides the default DOLFIN matrix class,
        based on the default DOLFIN linear algebra backend.

    .. cpp:function:: Matrix() : matrix(0)
    
        Create empty matrix

    .. cpp:function:: Matrix(const Matrix& A) : matrix(A.matrix->copy())
    
        Copy constructor

    .. cpp:function:: Matrix(uint M, uint N) : matrix(0)
    
        Create M x N matrix

    .. cpp:function:: const Matrix& operator= (const Matrix& A)
    
        Assignment operator

    .. cpp:function:: virtual GenericMatrix* instance()
    
        Return concrete instance / unwrap (non-const version)

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual Matrix* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator

    .. cpp:function:: virtual const GenericMatrix* instance() const
    
        Return concrete instance / unwrap (const version)

    .. cpp:function:: virtual const Matrix& operator*= (double a)
    
        Multiply matrix by given number

    .. cpp:function:: virtual const Matrix& operator/= (double a)
    
        Divide matrix by given number

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of matrix

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual std::tr1::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compressed storage data.
        See GenericMatrix for documentation.

    .. cpp:function:: virtual uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Add block of values

    .. cpp:function:: virtual void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericMatrix& A, bool same_nonzero_pattern)
    
        Add multiple of given matrix (AXPY operation)

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows, uint n, const uint* cols) const
    
        Get block of values

    .. cpp:function:: virtual void getrow(uint row, std::vector<uint>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row

    .. cpp:function:: virtual void ident(uint m, const uint* rows)
    
        Set given rows to identity matrix

    .. cpp:function:: virtual void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: virtual void resize(uint M, uint N)
    
        Resize matrix to M x N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Set block of values

    .. cpp:function:: virtual void setrow(uint row, const std::vector<uint>& columns, const std::vector<double>& values)
    
        Set values for given row

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual void zero(uint m, const uint* rows)
    
        Set given rows to zero

    .. cpp:function:: virtual ~Matrix()
    
        Destructor

