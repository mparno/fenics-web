.. Documentation for the header file dolfin/la/uBLASMatrix.h

.. _programmers_reference_cpp_la_ublasmatrix:

uBLASMatrix.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: uBLASMatrix

    *Parent class*
    
        * :cpp:class:`GenericMatrix`
        
        This class provides a simple matrix class based on uBLAS.
        It is a simple wrapper for a uBLAS matrix implementing the
        GenericMatrix interface.
        
        The interface is intentionally simple. For advanced usage,
        access the underlying uBLAS matrix and use the standard
        uBLAS interface which is documented at
        http://www.boost.org/libs/numeric/ublas/doc/index.htm.
        
        Developer note: specialised member functions must be
        inlined to avoid link errors.

    .. cpp:function:: Mat& mat()
    
        Return reference to uBLAS matrix (non-const version)

    .. cpp:function:: const Mat& mat() const
    
        Return reference to uBLAS matrix (const version)

    .. cpp:function:: const uBLASMatrix<Mat>& operator= (const uBLASMatrix<Mat>& A)
    
        Assignment operator

    .. cpp:function:: double operator() (uint i, uint j) const
    
        Access value of given entry

    .. cpp:function:: template <class E>
                                         explicit uBLASMatrix(const ublas::matrix_expression<E>& A) : Mat(A)
    
        Create matrix from given uBLAS matrix expression

    .. cpp:function:: template<class B>
                                        void solveInPlace(B& X)
    
        General uBLAS LU solver which accepts both vector and matrix right-hand sides

    .. cpp:function:: template<class Mat>
                       uBLASMatrix()
    
        Create empty matrix

    .. cpp:function:: uBLASMatrix(const uBLASMatrix& A)
    
        Copy constructor

    .. cpp:function:: uBLASMatrix(uint M, uint N)
    
        Create M x N matrix

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator

    .. cpp:function:: virtual const uBLASMatrix<Mat>& operator*= (double a)
    
        Multiply matrix by given number

    .. cpp:function:: virtual const uBLASMatrix<Mat>& operator/= (double a)
    
        Divide matrix by given number

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of matrix

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual std::tr1::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compresssed storage data
        See GenericMatrix for documentation.

    .. cpp:function:: virtual uBLASMatrix<Mat>* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Add block of values

    .. cpp:function:: virtual void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericMatrix& A,
                                        bool same_nonzero_pattern)
    
        Add multiple of given matrix (AXPY operation)

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows, uint n, const uint* cols) const
    
        Get block of values

    .. cpp:function:: virtual void getrow(uint row, std::vector<uint>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row

    .. cpp:function:: virtual void ident(uint m, const uint* rows)
    
        Set given rows to identity matrix

    .. cpp:function:: virtual void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: virtual void mult(const GenericVector& x, GenericVector& y) const
    
        Matrix-vector product, y = Ax

    .. cpp:function:: virtual void resize(uint M, uint N)
    
        Resize matrix to M x N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Set block of values

    .. cpp:function:: virtual void setrow(uint row_idx, const std::vector<uint>& columns, const std::vector<double>& values)
    
        Set values for given row

    .. cpp:function:: virtual void transpmult(const GenericVector& x, GenericVector& y) const
    
        Matrix-vector product, y = A^T x

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual void zero(uint m, const uint* rows)
    
        Set given rows to zero

    .. cpp:function:: virtual ~uBLASMatrix()
    
        Destructor

    .. cpp:function:: void compress()
    
        Compress matrix (eliminate all non-zeros from a sparse matrix)

    .. cpp:function:: void invert()
    
        Compute inverse of matrix

    .. cpp:function:: void lump(uBLASVector& m) const
    
        Lump matrix into vector m

    .. cpp:function:: void solve(uBLASVector& x, const uBLASVector& b) const
    
        Solve Ax = b out-of-place using uBLAS (A is not destroyed)

    .. cpp:function:: void solveInPlace(uBLASVector& x, const uBLASVector& b)
    
        Solve Ax = b in-place using uBLAS(A is destroyed)

