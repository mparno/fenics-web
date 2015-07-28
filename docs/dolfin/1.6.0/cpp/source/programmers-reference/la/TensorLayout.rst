
.. Documentation for the header file dolfin/la/TensorLayout.h

.. _programmers_reference_cpp_la_tensorlayout:

TensorLayout.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TensorLayout

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class described the size and possibly the sparsity of a
    (sparse) tensor. It is used by the linear algebra backends to
    initialise tensors.


    .. cpp:function:: TensorLayout(std::size_t primary_dim, bool sparsity_pattern)
    
        Create empty tensor layout


    .. cpp:function:: TensorLayout(const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, std::size_t primary_dim, std::size_t block_size, const std::vector<std::pair<std::size_t, std::size_t> >& ownership_range, bool sparsity_pattern)
    
        Create a tensor layout


    .. cpp:function:: void init(const MPI_Comm mpi_comm, const std::vector<std::size_t>& dims, std::size_t block_size, const std::vector<std::pair<std::size_t, std::size_t> >& ownership_range)
    
        Initialize tensor layout


    .. cpp:function:: std::size_t rank() const
    
        Return rank


    .. cpp:function:: std::size_t size(std::size_t i) const
    
        Return global size for dimension i (size of tensor, includes
        non-zeroes)


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local range for dimension dim


    .. cpp:function:: std::shared_ptr<GenericSparsityPattern> sparsity_pattern()
    
        Return sparsity pattern (possibly null)


    .. cpp:function:: std::shared_ptr<const GenericSparsityPattern> sparsity_pattern() const
    
        Return sparsity pattern (possibly null), const version


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


