.. Documentation for the header file dolfin/la/MTL4Matrix.h

.. _programmers_reference_cpp_la_mtl4matrix:

MTL4Matrix.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: implements a minimal backend for MTL4.

.. cpp:class:: MTL4Matrix

    *Parent class*
    
        * :cpp:class:`GenericMatrix`
        
    .. cpp:function:: MTL4Matrix()
    
        Create empty matrix

    .. cpp:function:: MTL4Matrix(const MTL4Matrix& A)
    
        Copy constuctor

    .. cpp:function:: MTL4Matrix(uint M, uint N)
    
        Create M x N matrix

    .. cpp:function:: MTL4Matrix(uint M, uint N, uint nz)
    
        Create M x N matrix with estimate of nonzeroes per row

    .. cpp:function:: const MTL4Matrix& operator= (const MTL4Matrix& A)
    
        Assignment operator

    .. cpp:function:: const mtl4_sparse_matrix& mat() const
    
        Return mtl4_sparse_matrix reference

    .. cpp:function:: virtual MTL4Matrix* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator

    .. cpp:function:: virtual const MTL4Matrix& operator*= (double a)
    
        Multiply matrix by given number

    .. cpp:function:: virtual const MTL4Matrix& operator/= (double a)
    
        Divide matrix by given number

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of matrix

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual std::tr1::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compresssed storage data
        See GenericMatrix for documentation.

    .. cpp:function:: virtual uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Add block of values

    .. cpp:function:: virtual void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericMatrix& A,  bool same_nonzero_pattern)
    
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

    .. cpp:function:: virtual ~MTL4Matrix()
    
        Destructor

