
.. Documentation for the header file dolfin/la/TrilinosPreconditioner.h

.. _programmers_reference_cpp_la_trilinospreconditioner:

TrilinosPreconditioner.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TrilinosPreconditioner

    *Parent class(es)*
    
        * :cpp:class:`GenericPreconditioner`
        
        * :cpp:class:`Variable`
        
    This class is a wrapper for configuring Epetra
    preconditioners. It does not own a preconditioner. It can take a
    EpetraKrylovSolver and set the preconditioner type and
    parameters.


    .. cpp:function:: explicit TrilinosPreconditioner(std::string method="default")
    
        Create Krylov solver for a particular method and preconditioner


    .. cpp:function:: void set(BelosLinearProblem& problem, const EpetraMatrix& P )
    
        Set the precondtioner and matrix used in preconditioner


    .. cpp:function:: void set_parameters(std::shared_ptr<const Teuchos::ParameterList> list)
    
        Set the Trilonos preconditioner parameters list


    .. cpp:function:: void set_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
    
        Set the Trilonos preconditioner parameters list (for use from
        Python)


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& null_space)
    
        Set basis for the null space of the operator. Setting this is
        critical to the performance of some preconditioners, e.g. ML.
        The vectors spanning the null space are copied.


    .. cpp:function:: std::string name() const
    
        Return preconditioner name


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void set_ml(BelosLinearProblem& problem, const Epetra_RowMatrix& P)
    
        Setup the ML precondtioner


