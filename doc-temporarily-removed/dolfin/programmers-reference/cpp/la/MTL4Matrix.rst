.. Documentation for the header file dolfin/la/MTL4Matrix.h

.. _programmers_reference_cpp_la_mtl4matrix:

MTL4Matrix.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MTL4Matrix

    *Parent class*
    
        * :cpp:class:`GenericMatrix`
        
    .. cpp:function:: MTL4Matrix()
    
        Create empty matrix

    .. cpp:function:: MTL4Matrix(const MTL4Matrix& A)
    
        Copy constuctor

    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: MTL4Matrix* copy() const
    
        Return copy of tensor

    .. cpp:function:: uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: std::pair<uint, uint> local_range(uint dim) const
    
        Return local ownership range

    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: void resize(uint M, uint N)
    
        Resize matrix to M x N

    .. cpp:function:: void resize(GenericVector& y, uint dim) const
    
        Resize vector y such that is it compatible with matrix for
        multuplication Ax = b (dim = 0 -> b, dim = 1 -> x) In parallel
        case, size and layout are important.

    .. cpp:function:: void get(double* block, uint m, const uint* rows, uint n, const uint* cols) const
    
        Get block of values

    .. cpp:function:: void set(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Set block of values

    .. cpp:function:: void add(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Add block of values

    .. cpp:function:: void axpy(double a, const GenericMatrix& A,  bool same_nonzero_pattern)
    
        Add multiple of given matrix (AXPY operation)

    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of matrix

    .. cpp:function:: void getrow(uint row, std::vector<uint>& columns, std::vector<double>& values) const
    
        Get non-zero values of given row

    .. cpp:function:: void setrow(uint row, const std::vector<uint>& columns, const std::vector<double>& values)
    
        Set values for given row

    .. cpp:function:: void zero(uint m, const uint* rows)
    
        Set given rows to zero

    .. cpp:function:: void ident(uint m, const uint* rows)
    
        Set given rows to identity matrix

    .. cpp:function:: const MTL4Matrix& operator*= (double a)
    
        Multiply matrix by given number

    .. cpp:function:: const MTL4Matrix& operator/= (double a)
    
        Divide matrix by given number

    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator

    .. cpp:function:: std::tr1::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compresssed storage data
        See GenericMatrix for documentation.

    .. cpp:function:: MTL4Matrix(uint M, uint N, uint nz)
    
        Create M x N matrix with estimate of nonzeroes per row

    .. cpp:function:: const mtl4_sparse_matrix& mat() const
    
        Return mtl4_sparse_matrix reference

    .. cpp:function:: const MTL4Matrix& operator= (const MTL4Matrix& A)
    
        Assignment operator

