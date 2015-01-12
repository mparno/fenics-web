
.. Documentation for the header file dolfin/adaptivity/GoalFunctional.h

.. _programmers_reference_cpp_adaptivity_goalfunctional:

GoalFunctional.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GoalFunctional

    *Parent class(es)*
    
        * :cpp:class:`Form`
        
    A :cpp:class:`GoalFunctional` is a :cpp:class:`Form` of rank 0 with an associated
    :cpp:class:`ErrorControl`.


    .. cpp:function:: GoalFunctional(std::size_t rank, std::size_t num_coefficients)
    
        Create :cpp:class:`GoalFunctional`
        
        *Arguments*
            rank (int)
                the rank of the functional (should be 0)
            num_coefficients (int)
                the number of coefficients in functional


    .. cpp:function:: void update_ec(const Form& a, const Form& L) = 0
    
        Update error control instance with given forms
        
        *Arguments*
            a (:cpp:class:`Form`)
                a bilinear form
            L (:cpp:class:`Form`)
                a linear form


