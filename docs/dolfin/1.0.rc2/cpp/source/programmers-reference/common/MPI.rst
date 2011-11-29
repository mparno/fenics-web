
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


.. cpp:class:: MPI

    This class provides utility functions for easy communcation with MPI.


    .. cpp:function:: static uint process_number()
    
        Return proccess number


    .. cpp:function:: static uint num_processes()
    
        Return number of processes


    .. cpp:function:: static bool is_broadcaster()
    
        Determine whether we should broadcast (based on current parallel policy)


    .. cpp:function:: static bool is_receiver()
    
        Determine whether we should receive (based on current parallel policy)


    .. cpp:function:: static void barrier()
    
        Set a barrier (synchronization point)


    .. cpp:function:: static void distribute(const std::vector<T>& in_values, const std::vector<uint>& destinations, std::vector<T>& out_values, std::vector<uint>& sources)
    
        Distribute local arrays on all processors according to given partition


    .. cpp:function:: static void distribute(const std::vector<T>& in_values, const std::vector<uint>& destinations, std::vector<T>& out_values)
    
        Distribute local arrays on all processors according to given partition


    .. cpp:function:: static void scatter(const std::vector<T>& in_values, T& out_value, uint sending_process=0)
    
        Scatter in_values[i] to process i


    .. cpp:function:: static uint global_offset(uint range, bool exclusive)
    
        Find global offset (index) (wrapper for MPI_(Ex)Scan with MPI_SUM as
        reduction op)


    .. cpp:function:: static void send_recv(const T& send_value, uint dest, T& recv_value, uint source)
    
        Send-receive and data


    .. cpp:function:: static std::pair<uint, uint> local_range(uint N)
    
        Return local range for local process, splitting [0, N - 1] into
        num_processes() portions of almost equal size


    .. cpp:function:: static std::pair<uint, uint> local_range(uint process, uint N)
    
        Return local range for given process, splitting [0, N - 1] into
        num_processes() portions of almost equal size


    .. cpp:function:: static std::pair<uint, uint> local_range(uint process, uint N, uint num_processes)
    
        Return local range for given process, splitting [0, N - 1] into
        num_processes portions of almost equal size


    .. cpp:function:: static uint index_owner(uint index, uint N)
    
        Return which process owns index (inverse of local_range)


