.. Documentation for the header file dolfin/la/Scalar.h

.. _programmers_reference_cpp_la_scalar:

Scalar.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Scalar

    *Parent class*
    
        * :cpp:class:`GenericTensor`
        
    This class represents a real-valued scalar quantity and
    implements the GenericTensor interface for scalars.

    .. cpp:function:: Scalar()
    
        Create zero scalar

    .. cpp:function:: void resize(uint rank, const uint* dims)
    
        Resize tensor to given dimensions

    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern)
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: Scalar* copy() const
    
        Return copy of tensor

    .. cpp:function:: uint rank() const
    
        Return tensor rank (number of dimensions)

    .. cpp:function:: uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: void get(double* block, const uint* num_rows, const uint * const * rows) const
    
        Get block of values

    .. cpp:function:: void set(const double* block, const uint* num_rows, const uint * const * rows)
    
        Set block of values

    .. cpp:function:: void add(const double* block, const uint* num_rows, const uint * const * rows)
    
        Add block of values

    .. cpp:function:: void add(const double* block, const std::vector<const std::vector<uint>* >& rows)
    
        Add block of values

    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: operator double() const
    
        Cast to real

    .. cpp:function:: const Scalar& operator= (double value)
    
        Assignment from real

    .. cpp:function:: LinearAlgebraFactory& factory() const
    
        Return a factory for the default linear algebra backend

    .. cpp:function:: double getval() const
    
        Get value

