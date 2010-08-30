.. Documentation for the header file dolfin/la/GenericMatrix.h

.. _programmers_reference_cpp_la_genericmatrix:

GenericMatrix.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GenericMatrix

    *Parent class*
    
        * :cpp:class:`GenericTensor`
        
    This class defines a common interface for matrices.

    .. cpp:function:: GenericMatrix* copy() const = 0
    
        Return copy of tensor

    .. cpp:function:: const GenericMatrix& operator*= (double a) = 0
    
        Multiply matrix by given number

    .. cpp:function:: const GenericMatrix& operator+= (const GenericMatrix& A)
    
        Add given matrix

    .. cpp:function:: const GenericMatrix& operator-= (const GenericMatrix& A)
    
        Subtract given matrix

    .. cpp:function:: const GenericMatrix& operator/= (double a) = 0
    
        Divide matrix by given number

    .. cpp:function:: const GenericMatrix& operator= (const GenericMatrix& x) = 0
    
        Assignment operator

    .. cpp:function:: double getitem(std::pair<uint, uint> ij) const
    
        Get value of given entry

    .. cpp:function:: double norm(std::string norm_type) const = 0
    
        Return norm of matrix

    .. cpp:function:: double operator() (uint i, uint j) const
    
        Get value of given entry

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: std::tr1::tuple<const std::size_t*, const std::size_t*, const double*, int> data() const
    
        Return pointers to underlying compresssed row/column storage data
        For compressed row storage, data = (row_pointer[#rows +1],
        column_index[#nz], matrix_values[#nz], nz)

    .. cpp:function:: uint rank() const
    
        Return tensor rank (number of dimensions)

    .. cpp:function:: uint size(uint dim) const = 0
    
        Return size of given dimension

    .. cpp:function:: void add(const double* block, const uint* num_rows, const uint * const * rows)
    
        Add block of values

    .. cpp:function:: void add(const double* block, uint m, const uint* rows, uint n, const uint* cols) = 0
    
        Add block of values

    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor

    .. cpp:function:: void axpy(double a, const GenericMatrix& A, bool same_nonzero_pattern) = 0
    
        Add multiple of given matrix (AXPY operation)

    .. cpp:function:: void get(double* block, const uint* num_rows, const uint * const * rows) const
    
        Get block of values

    .. cpp:function:: void get(double* block, uint m, const uint* rows, uint n, const uint* cols) const = 0
    
        Get block of values

    .. cpp:function:: void getrow(uint row, std::vector<uint>& columns, std::vector<double>& values) const = 0
    
        Get non-zero values of given row on local process

    .. cpp:function:: void ident(uint m, const uint* rows) = 0
    
        Set given rows to identity matrix

    .. cpp:function:: void ident_zeros()
    
        Insert one on the diagonal for all zero rows

    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern) = 0
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const = 0
    
        Matrix-vector product, y = Ax

    .. cpp:function:: void resize(uint M, uint N) = 0
    
        Resize matrix to  M x N

    .. cpp:function:: void resize(uint rank, const uint* dims)
    
        Resize tensor with given dimensions

    .. cpp:function:: void set(const double* block, const uint* num_rows, const uint * const * rows)
    
        Set block of values

    .. cpp:function:: void set(const double* block, uint m, const uint* rows, uint n, const uint* cols) = 0
    
        Set block of values

    .. cpp:function:: void setitem(std::pair<uint, uint> ij, double value)
    
        Set given entry to value. apply("insert") should be called before using
        using the object.

    .. cpp:function:: void setrow(uint row, const std::vector<uint>& columns, const std::vector<double>& values) = 0
    
        Set values for given row on local process

    .. cpp:function:: void transpmult(const GenericVector& x, GenericVector& y) const = 0
    
        Matrix-vector product, y = A^T x

    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void zero(uint m, const uint* rows) = 0
    
        Set given rows to zero

