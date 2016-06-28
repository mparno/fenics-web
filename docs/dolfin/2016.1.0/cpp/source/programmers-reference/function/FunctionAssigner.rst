
.. Documentation for the header file dolfin/function/FunctionAssigner.h

.. _programmers_reference_cpp_function_functionassigner:

FunctionAssigner.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FunctionAssigner

    This class facilitate assignments between Function and sub
    Functions. It builds and caches maps between compatible
    dofs. These maps are used in the assignment methods which
    perform the actual assignment. Optionally can a MeshFunction be
    passed together with a label, facilitating FunctionAssignment
    over sub domains.


    .. cpp:function:: FunctionAssigner(std::shared_ptr<const FunctionSpace> receiving_space, std::shared_ptr<const FunctionSpace> assigning_space)
    
        Create a FunctionAssigner between functions residing in the
        same type of FunctionSpace. One or both functions can be sub
        functions.
        
        *Arguments*
            receiving_space (:cpp:class:`FunctionSpace`)
                The function space of the receiving function
            assigning_space (:cpp:class:`FunctionSpace`)
                The function space of the assigning function


    .. cpp:function:: FunctionAssigner(std::vector<std::shared_ptr<const FunctionSpace> > receiving_spaces, std::shared_ptr<const FunctionSpace> assigning_space)
    
        Create a FunctionAssigner between one mixed function
        (assigning) and several functions (receiving). The number of
        receiving functions must sum up to the number of sub functions
        in the assigning mixed function. The sub spaces of the
        assigning mixed space must be of the same type ans size as the
        receiving spaces.
        
        *Arguments*
            receiving_spaces (std::vector<:cpp:class:`FunctionSpace`>)
                The receiving function spaces
            assigning_space (:cpp:class:`FunctionSpace`)
                The assigning function space


    .. cpp:function:: FunctionAssigner(std::shared_ptr<const FunctionSpace> receiving_space, std::vector<std::shared_ptr<const FunctionSpace> > assigning_spaces)
    
        Create a FunctionAssigner between several functions
        (assigning) and one mixed function (receiving). The number of
        sub functions in the assigning mixed function must sum up to
        the number of receiving functions. The sub spaces of the
        receiving mixed space must be of the same type ans size as the
        assigning spaces.
        
        *Arguments*
            receiving_space (std::shared_ptr<:cpp:class:`FunctionSpace`>)
                The receiving function space
            assigning_spaces (std::vector<std::shared_ptr<:cpp:class:`FunctionSpace`> >)
                The assigning function spaces


    .. cpp:function:: void assign(std::shared_ptr<Function> receiving_func, std::shared_ptr<const Function> assigning_func) const
    
        Assign one function to another
        
        *Arguments*
            receiving_func (std::shared_ptr<:cpp:class:`Function`>)
                The receiving function
            assigning_func (std::shared_ptr<:cpp:class:`Function`>)
                The assigning function


    .. cpp:function:: void assign(std::shared_ptr<Function> receiving_func, std::vector<std::shared_ptr<const Function> > assigning_funcs) const
    
        Assign several functions to sub functions of a mixed receiving
        function
        
        *Arguments*
            receiving_func (std::shared_ptr<:cpp:class:`Function`>)
                The receiving mixed function
            assigning_funcs (std::vector<std::shared_ptr<:cpp:class:`Function`> >)
                The assigning functions


    .. cpp:function:: void assign(std::vector<std::shared_ptr<Function> > receiving_funcs, std::shared_ptr<const Function> assigning_func) const
    
        Assign sub functions of a single mixed function to single
        receiving functions
        
        *Arguments*
            receiving_funcs (std::vector<std::shared_ptr<:cpp:class:`Function`> >)
                The receiving functions
            assigning_func (std::shared_ptr<:cpp:class:`Function`>)
                The assigning mixed function


    .. cpp:function:: std::size_t num_assigning_functions() const
    
        Return the number of assigning functions


    .. cpp:function:: std::size_t num_receiving_functions() const
    
        Return the number of receiving functions


