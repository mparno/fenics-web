
.. Documentation for the header file dolfin/fem/LinearVariationalProblem.h

.. _programmers_reference_cpp_fem_linearvariationalproblem:

LinearVariationalProblem.h
==========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LinearVariationalProblem

    *Parent class(es)*
    
        * :cpp:class:`Hierarchical<LinearVariationalProblem>`
        
    This class represents a linear variational problem:
    
    Find u in V such that
    
        a(u, v) = L(v)  for all v in V^,
    
    where V is the trial space and V^ is the test space.


    .. cpp:function:: LinearVariationalProblem(std::shared_ptr<const Form> a, std::shared_ptr<const Form> L, std::shared_ptr<Function> u, std::vector<std::shared_ptr<const DirichletBC>> bcs)
    
        Create linear variational problem with a list of boundary
        conditions


    .. cpp:function:: std::shared_ptr<const Form> bilinear_form() const
    
        Return bilinear form


    .. cpp:function:: std::shared_ptr<const Form> linear_form() const
    
        Return linear form


    .. cpp:function:: std::shared_ptr<Function> solution()
    
        Return solution variable


    .. cpp:function:: std::shared_ptr<const Function> solution() const
    
        Return solution variable (const version)


    .. cpp:function:: std::vector<std::shared_ptr<const DirichletBC>> bcs() const
    
        Return boundary conditions


    .. cpp:function:: std::shared_ptr<const FunctionSpace> trial_space() const
    
        Return trial space


    .. cpp:function:: std::shared_ptr<const FunctionSpace> test_space() const
    
        Return test space


