
.. Documentation for the header file dolfin/la/SparsityPattern.h

.. _programmers_reference_cpp_la_sparsitypattern:

SparsityPattern.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SparsityPattern

    This class implements a sparsity pattern data structure.  It is
    used by most linear algebra backends.


    .. cpp:function:: SparsityPattern(std::size_t primary_dim)
    
        Create empty sparsity pattern


    .. cpp:function:: SparsityPattern(MPI_Comm mpi_comm, std::vector<std::shared_ptr<const IndexMap>> index_maps, std::size_t primary_dim)
    
        Create sparsity pattern for a generic tensor


    .. cpp:function:: void init(MPI_Comm mpi_comm, std::vector<std::shared_ptr<const IndexMap>> index_maps)
    
        Initialize sparsity pattern for a generic tensor


    .. cpp:function:: void insert_global(dolfin::la_index i, dolfin::la_index j)
    
        Insert a global entry - will be fixed by apply()


    .. cpp:function:: void insert_global(const std::vector< ArrayView<const dolfin::la_index>>& entries)
    
        Insert non-zero entries using global indices


    .. cpp:function:: void insert_local(const std::vector< ArrayView<const dolfin::la_index>>& entries)
    
        Insert non-zero entries using local (process-wise) indices


    .. cpp:function:: void insert_full_rows_local(const std::vector<std::size_t>& rows)
    
        Insert full rows (or columns, according to primary dimension)
        using local (process-wise) indices. This must be called before
        any other sparse insertion occurs to avoid quadratic complexity
        of dense rows insertion


    .. cpp:function:: std::size_t rank() const
    
        Return rank


    .. cpp:function:: std::size_t primary_dim() const
    
        Return primary dimension (e.g., 0=row partition, 1=column
        partition)


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
        number of nonzeros per local row for off-diagonal block. If
        there is no off-diagonal pattern, the vector is resized to
        zero-length


    .. cpp:function:: void num_local_nonzeros(std::vector<std::size_t>& num_nonzeros) const
    
        Fill vector with number of nonzeros in local_range for
        dimension 0


    .. cpp:function:: void apply()
    
        Finalize sparsity pattern


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::vector<std::vector<std::size_t>> diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (diagonal). Options are
        'sorted' and 'unsorted'.


    .. cpp:function:: std::vector<std::vector<std::size_t>> off_diagonal_pattern(Type type) const
    
        Return underlying sparsity pattern (off-diagonal). Options are
        'sorted' and 'unsorted'. Empty vector is returned if there is no
        off-diagonal contribution.


