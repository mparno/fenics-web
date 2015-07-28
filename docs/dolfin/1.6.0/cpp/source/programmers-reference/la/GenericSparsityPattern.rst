
.. Documentation for the header file dolfin/la/GenericSparsityPattern.h

.. _programmers_reference_cpp_la_genericsparsitypattern:

GenericSparsityPattern.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericSparsityPattern

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    Base class (interface) for generic tensor sparsity patterns.
    Currently, this interface is mostly limited to matrices.


    .. cpp:function:: GenericSparsityPattern(std::size_t primary_dim)
    
        Create empty sparsity pattern


    .. cpp:function:: void init(const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, const std::vector<std::pair<std::size_t, std::size_t> >& local_range, const std::vector<ArrayView<const std::size_t> >& local_to_global, const std::vector<ArrayView<const int> >& off_process_owner, const std::vector<std::size_t>& block_sizes) = 0
    
        Initialize sparsity pattern for a generic tensor


    .. cpp:function:: void insert_global(const std::vector< ArrayView<const dolfin::la_index> >& entries) = 0
    
        Insert non-zero entries using global indices


    .. cpp:function:: void insert_local(const std::vector< ArrayView<const dolfin::la_index> >& entries) = 0
    
        Insert non-zero entries using local (process-wise) entries


    .. cpp:function:: std::size_t rank() const = 0
    
        Return rank


    .. cpp:function:: std::size_t primary_dim() const
    
        Return primary dimension (e.g., 0=row partition, 1=column
        partition)


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const = 0
    
        Return local range for dimension dim


    .. cpp:function:: std::size_t num_nonzeros() const = 0
    
        Return total number of nonzeros in local_range


    .. cpp:function:: void num_nonzeros_diagonal(std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros for diagonal block in
        local_range for primary dimension


    .. cpp:function:: void num_nonzeros_off_diagonal( std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros for off-diagonal block in
        local_range for primary dimension


    .. cpp:function:: void num_local_nonzeros(std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros in local_range for primary
        dimension


    .. cpp:function:: std::vector<std::vector<std::size_t> > diagonal_pattern(Type type) const = 0
    
        Return underlying sparsity pattern (diagonal). Options are
        'sorted' and 'unsorted'.


    .. cpp:function:: std::vector<std::vector<std::size_t> > off_diagonal_pattern(Type type) const = 0
    
        Return underlying sparsity pattern (off-diagonal). Options
        are 'sorted' and 'unsorted'.


    .. cpp:function:: void apply() = 0
    
        Finalize sparsity pattern


    .. cpp:function:: MPI_Comm mpi_comm() const = 0
    
        Return MPI communicator


