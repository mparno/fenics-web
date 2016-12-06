
.. Documentation for the header file dolfin/common/SubSystemsManager.h

.. _programmers_reference_cpp_common_subsystemsmanager:

SubSystemsManager.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SubSystemsManager

    This is a singleton class which manages the initialisation and
    finalisation of various sub systems, such as MPI and PETSc.


    .. cpp:function:: static SubSystemsManager& singleton()
    
        Singleton instance. Calling this ensures singleton instance of
        SubSystemsManager is initialized according to the "Construct
        on First Use" idiom.


    .. cpp:function:: static void init_mpi()
    
        Initialise MPI


    .. cpp:function:: static int init_mpi(int argc, char* argv[], int required_thread_level)
    
        Initialise MPI with required level of thread support


    .. cpp:function:: static void init_petsc()
    
        Initialize PETSc without command-line arguments


    .. cpp:function:: static void init_petsc(int argc, char* argv[])
    
        Initialize PETSc with command-line arguments. Note that PETSc
        command-line arguments may also be filtered and sent to PETSc
        by parameters.parse(argc, argv).


    .. cpp:function:: static void finalize()
    
        Finalize subsystems. This will be called by the destructor, but in
        special cases it may be necessary to call finalize() explicitly.


    .. cpp:function:: static bool responsible_mpi()
    
        Return true if DOLFIN initialised MPI (and is therefore responsible
        for finalization)


    .. cpp:function:: static bool responsible_petsc()
    
        Return true if DOLFIN initialised PETSc (and is therefore
        responsible for finalization)


    .. cpp:function:: static bool mpi_initialized()
    
        Check if MPI has been initialised (returns true if MPI has been
        initialised, even if it is later finalised)


    .. cpp:function:: static bool mpi_finalized()
    
        Check if MPI has been finalized (returns true if MPI has been
        finalised)


    .. cpp:function:: static PetscErrorCode PetscDolfinErrorHandler( MPI_Comm comm, int line, const char *fun, const char *file, PetscErrorCode n, PetscErrorType p, const char *mess, void *ctx)
    
        PETSc error handler. Logs everything known to DOLFIN logging
        system (with level TRACE) and stores the error message into
        pests_err_msg member.


