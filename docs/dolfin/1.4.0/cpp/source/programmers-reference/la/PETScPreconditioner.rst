
.. Documentation for the header file dolfin/la/PETScPreconditioner.h

.. _programmers_reference_cpp_la_petscpreconditioner:

PETScPreconditioner.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScPreconditioner

    *Parent class(es)*
    
        * :cpp:class:`PETScObject`
        
    This class is a wrapper for configuring PETSc
    preconditioners. It does not own a preconditioner. It can take a
    PETScKrylovSolver and set the preconditioner type and
    parameters.


    .. cpp:function:: explicit PETScPreconditioner(std::string type = "default")
    
        Create a particular preconditioner object


    .. cpp:function:: void set(PETScKrylovSolver& solver)
    
        Set the precondtioner type and parameters


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& near_nullspace)
    
        Set the (near) null space of the preconditioner operator
        (matrix). This is required for certain preconditioner types,
        e.g. smoothed aggregation multigrid


    .. cpp:function:: MatNullSpace near_nullspace() const
    
        Return the PETSc null space


    .. cpp:function:: void set_coordinates(const std::vector<double>& x, std::size_t dim)
    
        Set the coordinates of the operator (matrix) rows and
        geometric dimension d. This is can be used by required for
        certain preconditioners, e.g. ML. The input for this function
        can be generated using GenericDofMap::tabulate_all_dofs.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Rerturn a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


