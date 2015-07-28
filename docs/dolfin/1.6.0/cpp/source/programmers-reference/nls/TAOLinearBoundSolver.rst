
.. Documentation for the header file dolfin/nls/TAOLinearBoundSolver.h

.. _programmers_reference_cpp_nls_taolinearboundsolver:

TAOLinearBoundSolver.h
======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TAOLinearBoundSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
        * :cpp:class:`PETScObject`
        
    This class provides a bound constrained solver for a
    linear variational inequality defined by a matrix A and a vector b.
    It solves the problem:
    
    Find :math:`x_l\leq x\leq x_u` such that
    :math:`(Ax-b)\cdot (y-x)\geq 0,\; \forall x_l\leq y\leq x_u`
    
    It is a wrapper for the TAO bound constrained solver.
    
    *Example*
       .. code-block:: python
    
          # Assemble the linear system
          A, b = assemble_system(a, L, bc)
          # Define the constraints
          constraint_u = Constant(1.)
          constraint_l = Constant(0.)
          u_min = interpolate(constraint_l, V)
          u_max = interpolate(constraint_u, V)
          # Define the function to store the solution
          usol=Function(V)
          # Create the TAOLinearBoundSolver
          solver=TAOLinearBoundSolver("tao_gpcg","gmres")
          # Set some parameters
          solver.parameters["monitor_convergence"]=True
          solver.parameters["report"]=True
          # Solve the problem
          solver.solve(A, usol.vector(), b , u_min.vector(), u_max.vector())
          info(solver.parameters,True)
    


    .. cpp:function:: TAOLinearBoundSolver(const std::string method = "default", const std::string ksp_type = "default", const std::string pc_type = "default")
    
        Create TAO bound constrained solver


    .. cpp:function:: std::size_t solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b, const GenericVector& xl, const GenericVector& xu)
    
        Solve the linear variational inequality defined by A and b
        with xl =< x <= xu


    .. cpp:function:: std::size_t solve(const PETScMatrix& A, PETScVector& x, const PETScVector& b, const PETScVector& xl, const PETScVector& xu)
    
        Solve the linear variational inequality defined by A and b
        with xl =< x <= xu


    .. cpp:function:: void set_ksp( const std::string ksp_type = "default")
    
        Set PETSC Krylov Solver (ksp) used by TAO


    .. cpp:function:: static std::map<std::string, std::string> methods()
    
        Return a list of available Tao solver methods


    .. cpp:function:: static std::map<std::string, std::string> krylov_solvers()
    
        Return a list of available krylov solvers


    .. cpp:function:: static std::map<std::string, std::string> preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


