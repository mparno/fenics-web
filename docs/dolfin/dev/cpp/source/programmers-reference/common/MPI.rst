
.. Documentation for the header file dolfin/common/MPI.h

.. _programmers_reference_cpp_common_mpi:

MPI.h
=====

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MPICommunicator

    .. cpp:function:: MPICommunicator()
    
        Create communicator (copy of MPI_COMM_WORLD)


    .. cpp:function:: MPI_Comm& operator*()
    
        Dereference operator


.. cpp:class:: MPIInfo

.. cpp:class:: MPINonblocking

    .. cpp:function:: void send_recv(const T& send_value, unsigned int dest, T& recv_value, unsigned int source)
    
        Non-blocking send and receive


    .. cpp:function:: void send_recv(const T& send_value, unsigned int dest_tag, unsigned int dest, T& recv_value, unsigned int source_tag, unsigned int source)
    
        Non-blocking send and receive with tag


    .. cpp:function:: void wait_all()
    
        Wait for all requests to finish


.. cpp:class:: MPI

    This class provides utility functions for easy communcation
    with MPI.


    .. cpp:function:: static unsigned int process_number()
    
        Return proccess number


    .. cpp:function:: static unsigned int num_processes()
    
        Return number of processes


    .. cpp:function:: static bool is_broadcaster()
    
        Determine whether we should broadcast (based on current parallel
        policy)


    .. cpp:function:: static bool is_receiver()
    
        Determine whether we should receive (based on current parallel
        policy)


    .. cpp:function:: static void barrier()
    
        Set a barrier (synchronization point)


    .. cpp:function:: static void all_to_all(const std::vector<std::vector<T> >& in_values, std::vector<std::vector<T> >& out_values)
    
        Send in_values[p0] to process p0 and receive values from process
        p1 in out_values[p1]


    .. cpp:function:: static void distribute(const std::set<S> group, const std::map<S, T>& in_values_per_dest, std::map<S, T>& out_values_per_src)
    
        Distribute local arrays on a group of processes (typically
        neighbours from GenericDofMap::neighbours()). It is important
        that each process' group includes exactly the processes that
        has it in their groups, otherwise it will deadlock.


    .. cpp:function:: static void broadcast(T& value, unsigned int broadcaster=0)
    
        Broadcast value from broadcaster process to all processes


    .. cpp:function:: static void scatter(const std::vector<T>& in_values, T& out_value, unsigned int sending_process=0)
    
        Scatter in_values[i] to process i


    .. cpp:function:: static void gather(const T& in_value, std::vector<T>& out_values, unsigned int receiving_process=0)
    
        Gather values on one process (wrapper for boost::mpi::gather)


    .. cpp:function:: static void all_gather(const T& in_value, std::vector<T>& out_values)
    
        Gather values, one from each process (wrapper for boost::mpi::all_gather)


    .. cpp:function:: static T max(const T& value)
    
        Return global max value


    .. cpp:function:: static T min(const T& value)
    
        Return global min value


    .. cpp:function:: static T sum(const T& value)
    
        Sum values and return sum


    .. cpp:function:: static T all_reduce(const T& value, X op)
    
        All reduce


    .. cpp:function:: static std::size_t global_offset(std::size_t range, bool exclusive)
    
        Find global offset (index) (wrapper for MPI_(Ex)Scan with MPI_SUM as
        reduction op)


    .. cpp:function:: static void send_recv(const T& send_value, unsigned int dest, T& recv_value, unsigned int source)
    
        Send-receive data. Note that if the number of posted send-receives
        may differ between processes, another interface (such as
        MPINonblocking::send_recv) must be used since duplicating the
        communicator requires participation from all processes.


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(std::size_t N)
    
        Return local range for local process, splitting [0, N - 1] into
        num_processes() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(unsigned int process, std::size_t N)
    
        Return local range for given process, splitting [0, N - 1] into
        num_processes() portions of almost equal size


    .. cpp:function:: static std::pair<std::size_t, std::size_t> local_range(unsigned int process, std::size_t N, unsigned int num_processes)
    
        Return local range for given process, splitting [0, N - 1] into
        num_processes portions of almost equal size


    .. cpp:function:: static unsigned int index_owner(std::size_t index, std::size_t N)
    
        Return which process owns index (inverse of local_range)


