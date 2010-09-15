.. Documentation for the header file dolfin/la/GenericTensor.h

.. _programmers_reference_cpp_la_generictensor:

GenericTensor.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GenericTensor

    *Parent class*
    
        * :cpp:class:`virtual Variable`
        
    This class defines a common interface for arbitrary rank tensors.

    .. cpp:function:: void resize(uint rank, const uint* dims) = 0
    
        Resize tensor with given dimensions

    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern) = 0
    
        Initialize zero tensor using sparsity pattern

    .. cpp:function:: GenericTensor* copy() const = 0
    
        Return copy of tensor

    .. cpp:function:: uint rank() const = 0
    
        Return tensor rank (number of dimensions)

    .. cpp:function:: uint size(uint dim) const = 0
    
        Return size of given dimension

    .. cpp:function:: void get(double* block, const uint* num_rows, const uint * const * rows) const = 0
    
        Get block of values

    .. cpp:function:: void set(const double* block, const uint* num_rows, const uint * const * rows) = 0
    
        Set block of values

    .. cpp:function:: void add(const double* block, const uint* num_rows, const uint * const * rows) = 0
    
        Add block of values

    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: LinearAlgebraFactory& factory() const = 0
    
        Return linear algebra backend factory

    .. cpp:function:: template<class T> const T& down_cast() const
    
        Cast a GenericTensor to its derived class (const version)

    .. cpp:function:: template<class T> T& down_cast()
    
        Cast a GenericTensor to its derived class (non-const version)

    .. cpp:function:: template<class T> bool has_type() const
    
        Check whether the GenericTensor instance matches a specific type

    .. cpp:function:: const GenericTensor* instance() const
    
        Return concrete instance / unwrap (const version)

    .. cpp:function:: GenericTensor* instance()
    
        Return concrete instance / unwrap (non-const version)

    .. cpp:function:: const GenericTensor& operator= (const GenericTensor& x)
    
        Assignment (must be overloaded by subclass)

