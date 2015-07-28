
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


    .. cpp:function:: static unsigned int rank(MPI_Comm comm)
    
        Return process rank for the communicator


    .. cpp:function:: static unsigned int size(MPI_Comm comm)
    
        Return size of the group (number of processes) associated with
        the communicator


    .. cpp:function:: static bool is_broadcaster(MPI_Comm comm)
    
        Determine whether we should broadcast (based on current
        parallel policy)


    .. cpp:function:: static bool is_receiver(MPI_Comm comm)
    
        Determine whether we should receive (based on current parallel
        policy)


    .. cpp:function:: static void barrier(MPI_Comm comm)
    
        Set a barrier (synchronization point)


    .. cpp:function:: static void all_to_all(MPI_Comm comm, std::vector<std::vector<T> >& in_values, std::vector<std::vector<T> >& out_values)
    
        Send in_values[p0] to process p0 and receive values from
        process p1 in out_values[p1]


    .. cpp:function:: static void broadcast(MPI_Comm comm, std::vector<T>& value, unsigned int broadcaster=0)
    
        Broadcast vector of value from broadcaster to all processes


    .. cpp:function:: static void broadcast(MPI_Comm comm, T& value, unsigned int broadcaster=0)
    
        Broadcast single primitive from broadcaster to all processes


    .. cpp:function:: static void scatter(MPI_Comm comm, const std::vector<std::vector<T> >& in_values, std::vector<T>& out_value, unsigned int sending_process=0)
    
        Scatter vector in_values[i] to process i


    .. cpp:function:: static void scatter(MPI_Comm comm, const std::vector<T>& in_values, T& out_value, unsigned int sending_process=0)
    
        Scatter primitive in_values[i] to process i


    .. cpp:function:: static void gather(MPI_Comm comm, const std::vector<T>& in_values, std::vector<T>& out_values, unsigned int receiving_process=0)
    
        Gather values on one process


    .. cpp:function:: static void gather(MPI_Comm comm, const std::string& in_values, std::vector<std::string>& out_values, unsigned int receiving_process=0)
    
        Gather strings on one process


    .. cpp:function:: static void all_gather(MPI_Comm comm, const std::vector<T>& in_values, std::vector<T>& out_values)
    
        Gather values from all processes. Same data count from each
        process (wrapper for MPI_Allgather)


    .. cpp:function:: static void all_gather(MPI_Comm comm, const std::vector<T>& in_values, std::vector<std::vector<T> >& out_values)
    
        Gather values from each process (variable count per process)


    .. cpp:function:: static void all_gather(MPI_Comm comm, const T in_value, std::vector<T>& out_values)
    
        Gather values, one primitive from each process (MPI_Allgather)


    .. cpp:function:: static T max(MPI_Comm comm, const T& value)
    
        Return global max value


    .. cpp:function:: static T min(MPI_Comm comm, const T& value)
    
        Return global min value


    .. cpp:function:: static T sum(MPI_Comm comm, const T& value)
    
        Sum values and return sum


    .. cpp:function:: static T avg(MPI_Comm comm, const T& value)
    
        Return average across comm; implemented only for T == Table


    .. cpp:function:: static T all_reduce(MPI_Comm comm, const T& value, X op)
    
        All reduce


    .. cpp:function:: static std::size_t global_offset(MPI_Comm comm, std::size_t range, bool exclusive)
    
        Find global offset (index) (wrapper for MPI_(Ex)Scan with
        MPI_SUM as reduction op)


    .. cpp:function:: static void send_recv(MPI_Comm comm, const std::vector<T>& send_value, unsigned int dest, int send_tag, std::vector<T>& recv_value, unsigned int source, int recv_tag)
    
        Send-receive data between processes (blocking)


    .. cpp:function:: static void send_recv(MPI_Comm comm, const std::vector<T>& send_value, unsigned int dest, std::vector<T>& recv_value, unsigned int source)
    
        Send-receive data between processes


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(MPI_Comm comm, std::size_t N)
    
        Return local range for local process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(MPI_Comm comm, unsigned int process, std::size_t N)
    
        Return local range for given process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> compute_local_range(unsigned int process, std::size_t N, unsigned int size)
    
        Return local range for given process, splitting [0, N - 1] into
        size() portions of almost equal size


    .. cpp:function:: static unsigned int index_owner(MPI_Comm comm, std::size_t index, std::size_t N)
    
        Return which process owns index (inverse of local_range)


    .. cpp:function:: static MPI_Op MPI_AVG()
    
        Return average reduction operation; recognized by
        all_reduce(MPI_Comm, Table&, MPI_Op)


