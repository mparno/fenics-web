
.. Documentation for the header file dolfin/fem/VariationalProblem.h

.. _programmers_reference_cpp_fem_variationalproblem:

VariationalProblem.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: VariationalProblem

    This class is deprecated and is only here to give an informative error
    message to users about the new interface.


    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1)
    
        Deprecated


    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1, const BoundaryCondition& bc)
    
        Deprecated


    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1, const std::vector<const BoundaryCondition*> bcs)
    
        Deprecated


    .. cpp:function:: VariationalProblem(boost::shared_ptr<const Form> form_0, boost::shared_ptr<const Form> form_1, std::vector<boost::shared_ptr<const BoundaryCondition> > bcs)
    
        Deprecated


    .. cpp:function:: void solve(Function& u) const
    
        Deprecated


    .. cpp:function:: void solve(Function& u0, Function& u1) const
    
        Deprecated


    .. cpp:function:: void solve(Function& u0, Function& u1, Function& u2) const
    
        Deprecated


    .. cpp:function:: void solve(Function& u, const double tol, GoalFunctional& M) const
    
        Deprecated


    .. cpp:function:: void solve(Function& u, const double tol, Form& M, ErrorControl& ec) const
    
        Deprecated


