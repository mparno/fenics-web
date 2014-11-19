
.. Documentation for the header file dolfin/function/assign.h

.. _programmers_reference_cpp_function_assign:

assign.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void assign(std::shared_ptr<Function> receiving_func, std::shared_ptr<const Function> assigning_func)

    Assign one function to another. The functions must reside in the
    same type of FunctionSpace. One or both functions can be sub
    functions.
    
    *Arguments*
        receiving_func (std::shared_ptr<:cpp:class:`Function`>)
            The recieving function
        assigning_func (std::shared_ptr<:cpp:class:`Function`>)
            The assigning function


.. cpp:function:: void assign(std::shared_ptr<Function> receiving_func, std::vector<std::shared_ptr<const Function> > assigning_funcs)

    Assign several functions to sub functions of a mixed receiving
    function. The number of receiving functions must sum up to the
    number of sub functions in the assigning mixed function. The sub
    spaces of the assigning mixed space must be of the same type ans
    size as the receiving spaces.


.. cpp:function:: void assign(std::vector<std::shared_ptr<Function> > receiving_funcs, std::shared_ptr<const Function> assigning_func)

    Assign sub functions of a single mixed function to single
    receiving functions. The number of sub functions in the
    assigning mixed function must sum up to the number of receiving
    functions. The sub spaces of the receiving mixed space must be
    of the same type ans size as the assigning spaces.


