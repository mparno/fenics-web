.. Documentation for the header file dolfin/la/SparsityPattern.h

.. _programmers_reference_cpp_la_sparsitypattern:

SparsityPattern.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: SparsityPattern

    *Parent class*
    
        * :cpp:class:`GenericSparsityPattern`
        
    This class implements the GenericSparsityPattern interface.
    It is used by most linear algebra backends.

    .. cpp:function:: SparsityPattern()
    
        Create empty sparsity pattern

    .. cpp:function:: void init(uint rank, const uint* dims)
    
        Initialize sparsity pattern for a generic tensor

    .. cpp:function:: void insert(const uint* num_rows, const uint * const * rows)
    
        Insert non-zero entries

    .. cpp:function:: uint rank() const
    
        Return rank

    .. cpp:function:: uint size(uint i) const
    
        Return global size for dimension i

    .. cpp:function:: std::pair<uint, uint> local_range(uint dim) const
    
        Return local range for dimension dim

    .. cpp:function:: uint num_nonzeros() const
    
        Return total number of nonzeros in local_range for dimension 0

    .. cpp:function:: void num_nonzeros_diagonal(uint* num_nonzeros) const
    
        Fill array with number of nonzeros for diagonal block in local_range for dimension 0
        For matrices, fill array with number of nonzeros per local row for diagonal block

    .. cpp:function:: void num_nonzeros_off_diagonal(uint* num_nonzeros) const
    
        Fill array with number of nonzeros for off-diagonal block in local_range for dimension 0
        For matrices, fill array with number of nonzeros per local row for off-diagonal block

    .. cpp:function:: void apply()
    
        Finalize sparsity pattern

    .. cpp:function:: std::string str() const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: std::vector<std::vector<uint> > diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (diagonal). Options are
        'sorted' and 'unsorted'.

    .. cpp:function:: std::vector<std::vector<uint> > off_diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (off-diagional). Options are
        'sorted' and 'unsorted'.

