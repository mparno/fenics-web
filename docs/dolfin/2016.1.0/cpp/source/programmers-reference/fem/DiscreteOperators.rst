
.. Documentation for the header file dolfin/fem/DiscreteOperators.h

.. _programmers_reference_cpp_fem_discreteoperators:

DiscreteOperators.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DiscreteOperators

    This class computes discrete gradient operators (matrices) that
    map derivatives of finite element functions into other finite
    element spaces. An example of where discrete gradient operators
    are required is the creation of algebraic multigrid solvers for
    H(curl) and H(div) problems.
    NOTE: This class is highly experimental and likely to change. It
    will eventually be expanded to provide the discrete curl and
    divergence.


    .. cpp:function:: static std::shared_ptr<GenericMatrix> build_gradient(const FunctionSpace& V0, const FunctionSpace& V1)
    
        Build the discrete gradient operator A that takes a w \in H^1
        (P1, nodal Lagrange) to v \in H(curl) (lowest order Nedelec),
        i.e. v = Aw. V0 is the H(curl) space, and V1 is the P1
        Lagrange space.


