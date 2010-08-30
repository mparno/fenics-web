.. Documentation for the header file dolfin/main/MPI.h

.. _programmers_reference_cpp_main_mpi:

MPI.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MPICommunicator

    .. cpp:function:: MPICommunicator()
    
        Create communicator (copy of MPI_COMM_WORLD)

    .. cpp:function:: MPI_Comm& operator*()
    
        Dereference operator

.. cpp:class:: MPI

    This class provides utility functions for easy communcation with MPI.

    .. cpp:function:: static bool is_broadcaster()
    
        Determine whether we should broadcast (based on current parallel policy)

    .. cpp:function:: static bool is_receiver()
    
        Determine whether we should receive (based on current parallel policy)

    .. cpp:function:: static double sum(double value)
    
        Sum values and return sum

    .. cpp:function:: static std::pair<uint, uint> local_range(uint N)
    
        Return local range for local process, splitting [0, N - 1] into
         num_processes() portions of almost equal size

    .. cpp:function:: static std::pair<uint, uint> local_range(uint process, uint N)
    
        Return local range for given process, splitting [0, N - 1] into
        num_processes() portions of almost equal size

    .. cpp:function:: static std::vector<uint> gather(uint value)
    
        Gather values, one from each process (wrapper for MPI_Allgather)

    .. cpp:function:: static uint global_maximum(uint size)
    
        Find global max value (wrapper for MPI_Allredue with MPI_MAX as
        reduction op)

    .. cpp:function:: static uint global_offset(uint range, bool exclusive)
    
        Find global offset (index) (wrapper for MPI_(Ex)Scan with MPI_SUM as
        reduction op)

    .. cpp:function:: static uint index_owner(uint index, uint N)
    
        Return which process owns index (inverse of local_range)

    .. cpp:function:: static uint num_processes()
    
        Return number of processes

    .. cpp:function:: static uint process_number()
    
        Return proccess number

    .. cpp:function:: static uint send_recv(double* send_buffer, uint send_size, uint dest,
                                            double* recv_buffer, uint recv_size, uint source)
    
        Send-receive and return number of received values (wrapper for MPI_Sendrecv)

    .. cpp:function:: static uint send_recv(uint* send_buffer, uint send_size, uint dest,
                                            uint* recv_buffer, uint recv_size, uint source)
    
        Send-receive and return number of received values (wrapper for MPI_Sendrecv)

    .. cpp:function:: static void barrier()
    
        Set a barrier (synchronization point)

    .. cpp:function:: static void distribute(std::vector<double>& values, std::vector<uint>& partition)
    
        Distribute local arrays on all processors according to given partition

    .. cpp:function:: static void distribute(std::vector<uint>& values, std::vector<uint>& partition)
    
        Distribute local arrays on all processors according to given partition

    .. cpp:function:: static void gather(std::vector<double>& values)
    
        Gather values, one from each process (wrapper for MPI_Allgather)

    .. cpp:function:: static void gather(std::vector<uint>& values)
    
        Gather values, one from each process (wrapper for MPI_Allgather)

    .. cpp:function:: static void scatter(std::vector<std::vector<double> >& values, uint sending_process=0)
    
        Scatter values (wrapper for MPI_Scatterv)

    .. cpp:function:: static void scatter(std::vector<std::vector<uint> >& values, uint sending_process=0)
    
        Scatter values (wrapper for MPI_Scatterv)

    .. cpp:function:: static void scatter(std::vector<uint>& values, uint sending_process=0)
    
        Scatter values, one to each process

    .. cpp:function:: template<class T> static void gather_all(const T& in_value,
                                                               std::vector<T>& out_values)
    
        Gather values, one from each process (wrapper for boost::mpi::all_gather)

