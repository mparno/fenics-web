.. Documentation for the header file dolfin/la/GenericSparsityPattern.h

.. _programmers_reference_cpp_la_genericsparsitypattern:

GenericSparsityPattern.h
========================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GenericSparsityPattern

    *Parent class*
    
        * :cpp:class:`Variable`
        
    Base class (interface) for generic tensor sparsity patterns.
    Currently, this interface is mostly limited to matrices.

    .. cpp:function:: GenericSparsityPattern()
    
        Create empty sparsity pattern

    .. cpp:function:: void init(uint rank, const uint* dims) = 0
    
        Initialize sparsity pattern for a generic tensor

    .. cpp:function:: void insert(const uint* num_rows, const uint * const * rows) = 0
    
        Insert non-zero entries

    .. cpp:function:: uint rank() const = 0
    
        Return rank

    .. cpp:function:: uint size(uint i) const = 0
    
        Return global size for dimension i

    .. cpp:function:: std::pair<uint, uint> local_range(uint dim) const = 0
    
        Return local range for dimension dim

    .. cpp:function:: uint num_nonzeros() const = 0
    
        Return total number of nonzeros in local_range for dimension 0

    .. cpp:function:: void num_nonzeros_diagonal(uint* num_nonzeros) const = 0
    
        Fill array with number of nonzeros for diagonal block in local_range for dimension 0

    .. cpp:function:: void num_nonzeros_off_diagonal(uint* num_nonzeros) const = 0
    
        Fill array with number of nonzeros for off-diagonal block in local_range for dimension 0

    .. cpp:function:: void apply() = 0
    
        Finalize sparsity pattern

