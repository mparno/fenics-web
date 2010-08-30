.. Documentation for the header file dolfin/la/STLMatrix.h

.. _programmers_reference_cpp_la_stlmatrix:

STLMatrix.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: STLMatrix

    *Parent class*
    
        * :cpp:class:`GenericMatrix`
        
    Simple implementation of a GenericMatrix for experimenting
    with new assembly. Not sure this will be used later but it
    might be useful.

    .. cpp:function:: STLMatrix()
    
        Create empty matrix

    .. cpp:function:: STLMatrix(const STLMatrix& A)
    
        Copy constructor

    .. cpp:function:: STLMatrix(uint M, uint N)
    
        Create M x N matrix

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        --- Specialized matrix functions ---
        Return linear algebra backend factory

    .. cpp:function:: virtual STLMatrix* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual const GenericMatrix& operator= (const GenericMatrix& A)
    
        Assignment operator

    .. cpp:function:: virtual const STLMatrix& operator*= (double a)
    
        Multiply matrix by given number

    .. cpp:function:: virtual const STLMatrix& operator/= (double a)
    
        Divide matrix by given number

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of matrix

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

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
    
        --- Implementation of the GenericTensor interface ---
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: virtual void resize(uint M, uint N)
    
        Initialize M x N matrix

    .. cpp:function:: virtual void resize(uint rank, const uint* dims, bool reset)
    
        Resize tensor of given rank and dimensions

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows, uint n, const uint* cols)
    
        Set block of values

    .. cpp:function:: virtual void setrow(uint row, const std::vector<uint>& columns, const std::vector<double>& values)
    
        Set values for given row

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual void zero(uint m, const uint* rows)
    
        Set given rows to zero

    .. cpp:function:: virtual ~STLMatrix()
    
        Destructor

