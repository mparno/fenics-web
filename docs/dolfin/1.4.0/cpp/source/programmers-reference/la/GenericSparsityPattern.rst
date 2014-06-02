
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


    .. cpp:function:: void init(const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, const std::vector<std::pair<std::size_t, std::size_t> >& local_range, const std::vector<const boost::unordered_map<std::size_t, unsigned int>* > off_process_owner) = 0
    
        Initialize sparsity pattern for a generic tensor


    .. cpp:function:: void insert(const std::vector<const std::vector<dolfin::la_index>* >& entries) = 0
    
        Insert non-zero entries


    .. cpp:function:: void add_edges(const std::pair<dolfin::la_index, std::size_t>& vertex, const std::vector<dolfin::la_index>& edges) = 0
    
        Add edges (vertex = [index, owning process])


    .. cpp:function:: std::size_t rank() const = 0
    
        Return rank


    .. cpp:function:: std::size_t primary_dim() const
    
        Return primary dimension (e.g., 0=row partition, 1=column partition)


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const = 0
    
        Return local range for dimension dim


    .. cpp:function:: std::size_t num_nonzeros() const = 0
    
        Return total number of nonzeros in local_range


    .. cpp:function:: void num_nonzeros_diagonal(std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros for diagonal block in
        local_range for primary dimemsion


    .. cpp:function:: void num_nonzeros_off_diagonal(std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros for off-diagonal block in
        local_range for primary dimemsion


    .. cpp:function:: void num_local_nonzeros(std::vector<std::size_t>& num_nonzeros) const = 0
    
        Fill vector with number of nonzeros in local_range for
        primary dimemsion


    .. cpp:function:: std::vector<std::vector<std::size_t> > diagonal_pattern(Type type) const = 0
    
        Return underlying sparsity pattern (diagonal). Options are
        'sorted' and 'unsorted'.


    .. cpp:function:: std::vector<std::vector<std::size_t> > off_diagonal_pattern(Type type) const = 0
    
        Return underlying sparsity pattern (off-diagional). Options are
        'sorted' and 'unsorted'.


    .. cpp:function:: void get_edges(std::size_t vertex, std::vector<dolfin::la_index>& edges) const = 0
    
        Fill vector with edges for given vertex


    .. cpp:function:: void apply() = 0
    
        Finalize sparsity pattern


    .. cpp:function:: MPI_Comm mpi_comm() const = 0
    
        Return MPI communicator


