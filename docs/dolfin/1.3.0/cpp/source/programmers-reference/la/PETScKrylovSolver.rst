
.. Documentation for the header file dolfin/la/PETScKrylovSolver.h

.. _programmers_reference_cpp_la_petsckrylovsolver:

PETScKrylovSolver.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
        * :cpp:class:`PETScObject`
        
    This class implements Krylov methods for linear systems
    of the form Ax = b. It is a wrapper for the Krylov solvers
    of PETSc.


    .. cpp:function:: PETScKrylovSolver(std::string method = "default", std::string preconditioner = "default")
    
        Create Krylov solver for a particular method and names
        preconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method, PETScPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method, boost::shared_ptr<PETScPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: PETScKrylovSolver(std::string method, PETScUserPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method, boost::shared_ptr<PETScUserPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: explicit PETScKrylovSolver(boost::shared_ptr<KSP> ksp)
    
        Create solver from given PETSc KSP pointer


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operator(const boost::shared_ptr<const PETScBaseMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericLinearOperator> A, const boost::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: void set_operators(const boost::shared_ptr<const PETScBaseMatrix> A, const boost::shared_ptr<const PETScBaseMatrix> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& nullspace)
    
        Set null space of the operator (matrix). This is used to solve
        singular systems


    .. cpp:function:: const PETScBaseMatrix& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const PETScBaseMatrix& A, PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<KSP> ksp() const
    
        Return PETSc KSP pointer


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: void set_options_prefix(std::string prefix)
    
        Set options prefix


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


