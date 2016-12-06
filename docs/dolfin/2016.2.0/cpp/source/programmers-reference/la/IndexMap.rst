
.. Documentation for the header file dolfin/la/IndexMap.h

.. _programmers_reference_cpp_la_indexmap:

IndexMap.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: IndexMap

    This class represents the distribution index arrays across
    processes. An index array is a contiguous collection of N+1
    indices [0, 1, . . ., N] that are distributed across processes M
    processes. On a given process, the IndexMap stores a portion of
    the index set using local indices [0, 1, . . . , n], and a map
    from the local indices to a unique global index.


.. cpp:class:: MapSize

    .. cpp:function:: IndexMap()
    
        Constructor


    .. cpp:function:: explicit IndexMap(MPI_Comm mpi_comm)
    
        Index map with no data


    .. cpp:function:: IndexMap(MPI_Comm mpi_comm, std::size_t local_size, std::size_t block_size)
    
        Index map with local size on each process. This constructor
        is collective


    .. cpp:function:: void init(std::size_t local_size, std::size_t block_size)
    
        Initialise with number of local entries and block size. This
        function is collective


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range() const
    
        Local range of indices


    .. cpp:function:: std::size_t size(MapSize type) const
    
        Get number of local indices of type MapSize::OWNED,
        MapSize::UNOWNED, MapSize::ALL or MapSize::GLOBAL


    .. cpp:function:: const std::vector<std::size_t>& local_to_global_unowned() const
    
        Get local to global map for unowned indices
        (local indexing beyond end of local range)


    .. cpp:function:: std::size_t local_to_global(std::size_t i) const
    
        Get global index of local index i


    .. cpp:function:: void set_local_to_global(const std::vector<std::size_t>& indices)
    
        Set local_to_global map for unowned indices (beyond end of local
        range). Computes and stores off-process owner array.


    .. cpp:function:: const std::vector<int>& off_process_owner() const
    
        Get off process owner for unowned indices


    .. cpp:function:: int block_size() const
    
        Get block size


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


