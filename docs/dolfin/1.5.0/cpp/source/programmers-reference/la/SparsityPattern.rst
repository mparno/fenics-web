
.. Documentation for the header file dolfin/la/SparsityPattern.h

.. _programmers_reference_cpp_la_sparsitypattern:

SparsityPattern.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SparsityPattern

    *Parent class(es)*
    
        * :cpp:class:`GenericSparsityPattern`
        
    This class implements the GenericSparsityPattern interface.  It
    is used by most linear algebra backends.


    .. cpp:function:: SparsityPattern(std::size_t primary_dim)
    
        Create empty sparsity pattern


    .. cpp:function:: SparsityPattern( const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, const std::vector<std::pair<std::size_t, std::size_t> >& ownership_range, const std::vector<const std::vector<std::size_t>* > local_to_global, const std::vector<const std::vector<int>* > off_process_owner, const std::vector<std::size_t>& block_sizes, const std::size_t primary_dim)
    
        Create sparsity pattern for a generic tensor


    .. cpp:function:: void init( const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, const std::vector<std::pair<std::size_t, std::size_t> >& ownership_range, const std::vector<const std::vector<std::size_t>* > local_to_global, const std::vector<const std::vector<int>* > off_process_owner, const std::vector<std::size_t>& block_sizes)
    
        Initialize sparsity pattern for a generic tensor


    .. cpp:function:: void insert_global(const std::vector< const std::vector<dolfin::la_index>* >& entries)
    
        Insert non-zero entries using global indices


    .. cpp:function:: void insert_local(const std::vector< const std::vector<dolfin::la_index>* >& entries)
    
        Insert non-zero entries using local (process-wise) indices


    .. cpp:function:: std::size_t rank() const
    
        Return rank


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local range for dimension dim


    .. cpp:function:: std::size_t num_nonzeros() const
    
        Return number of local nonzeros


    .. cpp:function:: void num_nonzeros_diagonal(std::vector<std::size_t>& num_nonzeros) const
    
        Fill array with number of nonzeros for diagonal block in
        local_range for dimension 0. For matrices, fill array with
        number of nonzeros per local row for diagonal block


    .. cpp:function:: void num_nonzeros_off_diagonal(std::vector<std::size_t>& num_nonzeros) const
    
        Fill array with number of nonzeros for off-diagonal block in
        local_range for dimension 0. For matrices, fill array with
        number of nonzeros per local row for off-diagonal block


    .. cpp:function:: void num_local_nonzeros(std::vector<std::size_t>& num_nonzeros) const
    
        Fill vector with number of nonzeros in local_range for
        dimension 0


    .. cpp:function:: void apply()
    
        Finalize sparsity pattern


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::vector<std::vector<std::size_t> > diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (diagonal). Options are
        'sorted' and 'unsorted'.


    .. cpp:function:: std::vector<std::vector<std::size_t> > off_diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (off-diagonal). Options are
        'sorted' and 'unsorted'.


