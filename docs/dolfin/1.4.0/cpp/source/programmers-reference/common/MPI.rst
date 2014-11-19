
.. Documentation for the header file dolfin/common/MPI.h

.. _programmers_reference_cpp_common_mpi:

MPI.h
=====

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MPIInfo

.. cpp:class:: MPI

    This class provides utility functions for easy communication
    with MPI and handles cases when DOLFIN is not configured with
    MPI.


    .. cpp:function:: static unsigned int process_number()
    
        Return process rank (uses MPI_COMM_WORLD)
        Warning: This function is deprecated. Use dolfin::MPI::rank


    .. cpp:function:: static unsigned int num_processes()
    
        Return number of processes for MPI_COMM_WORLD.
        Warning: This function is deprecated. Use dolfin::MPI::size.


    .. cpp:function:: static unsigned int rank(const MPI_Comm comm)
    
        Return process rank for the communicator


    .. cpp:function:: static unsigned int size(const MPI_Comm comm)
    
        Return size of the group (number of processes) associated with
        the communicator


    .. cpp:function:: static bool is_broadcaster(const MPI_Comm comm)
    
        Determine whether we should broadcast (based on current
        parallel policy)


    .. cpp:function:: static bool is_receiver(const MPI_Comm comm)
    
        Determine whether we should receive (based on current parallel
        policy)


    .. cpp:function:: static void barrier(const MPI_Comm comm)
    
        Set a barrier (synchronization point)


    .. cpp:function:: static void all_to_all(const MPI_Comm comm, std::vector<std::vector<T> >& in_values, std::vector<std::vector<T> >& out_values)
    
        Send in_values[p0] to process p0 and receive values from
        process p1 in out_values[p1]


    .. cpp:function:: static void broadcast(const MPI_Comm comm, std::vector<T>& value, unsigned int broadcaster=0)
    
        Broadcast vector of value from broadcaster to all processes


    .. cpp:function:: static void broadcast(const MPI_Comm comm, T& value, unsigned int broadcaster=0)
    
        Broadcast single primitive from broadcaster to all processes


    .. cpp:function:: static void scatter(const MPI_Comm comm, const std::vector<std::vector<T> >& in_values, std::vector<T>& out_value, unsigned int sending_process=0)
    
        Scatter vector in_values[i] to process i


    .. cpp:function:: static void scatter(const MPI_Comm comm, const std::vector<T>& in_values, T& out_value, unsigned int sending_process=0)
    
        Scatter primitive in_values[i] to process i


    .. cpp:function:: static void gather(const MPI_Comm comm, const std::vector<T>& in_values, std::vector<T>& out_values, unsigned int receiving_process=0)
    
        Gather values on one process


    .. cpp:function:: static void gather(const MPI_Comm comm, const std::string& in_values, std::vector<std::string>& out_values, unsigned int receiving_process=0)
    
        Gather strings on one process


    .. cpp:function:: static void all_gather(const MPI_Comm comm, const std::vector<T>& in_values, std::vector<T>& out_values)
    
        Gather values from all proceses. Same data count from each
        process (wrapper for MPI_Allgather)


    .. cpp:function:: static void all_gather(const MPI_Comm comm, const std::vector<T>& in_values, std::vector<std::vector<T> >& out_values)
    
        Gather values from each process (variable count per process)


    .. cpp:function:: static void all_gather(const MPI_Comm comm, const T& in_value, std::vector<T>& out_values)
    
        Gather values, one primitive from each process (MPI_Allgather)


    .. cpp:function:: static T max(const MPI_Comm comm, const T& value)
    
        Return global max value


    .. cpp:function:: static T min(const MPI_Comm comm, const T& value)
    
        Return global min value


    .. cpp:function:: static T sum(const MPI_Comm comm, const T& value)
    
        Sum values and return sum


    .. cpp:function:: static T all_reduce(const MPI_Comm comm, const T& value, X op)
    
        All reduce


    .. cpp:function:: static std::size_t global_offset(const MPI_Comm comm, std::size_t range, bool exclusive)
    
        Find global offset (index) (wrapper for MPI_(Ex)Scan with
        MPI_SUM as reduction op)


    .. cpp:function:: static void send_recv(const MPI_Comm comm, const std::vector<T>& send_value, unsigned int dest, int send_tag, std::vector<T>& recv_value, unsigned int source, int recv_tag)
    
        Send-receive data between processes (blocking)


    .. cpp:function:: static void send_recv(const MPI_Comm comm, const std::vector<T>& send_value, unsigned int dest, std::vector<T>& recv_value, unsigned int source)
    
        Send-receive data between processes


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(const MPI_Comm comm, std::size_t N)
    
        Return local range for local process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(const MPI_Comm comm, unsigned int process, std::size_t N)
    
        Return local range for given process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> compute_local_range(unsigned int process, std::size_t N, unsigned int size)
    
        Return local range for given process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static unsigned int index_owner(const MPI_Comm comm, std::size_t index, std::size_t N)
    
        Return which process owns index (inverse of local_range)


