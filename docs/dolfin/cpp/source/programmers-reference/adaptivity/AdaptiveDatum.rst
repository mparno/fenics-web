
.. Documentation for the header file dolfin/adaptivity/AdaptiveDatum.h

.. _programmers_reference_cpp_adaptivity_adaptivedatum:

AdaptiveDatum.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: AdaptiveDatum

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    An :cpp:class:`AdaptiveDatum` is a storage unit for data created in an
    adaptive process.


    .. cpp:function:: AdaptiveDatum(const uint refinement_level, const uint num_dofs, const uint num_cells, const double error_estimate, const double tolerance, const double functional_value)
    
        Create adaptive datum
        
        *Arguments*
            refinement_level (unsigned int)
                the number of refinements relative to coarset mesh
            num_dofs (unsigned int)
                dimension of discrete solution space
            num_cells (unsigned int)
                number of cells in mesh
            error_estimate (double)
                error estimate
            tolerance (double)
                error (or num_dofs) tolerance


    .. cpp:function:: void store(std::string filename) const
    
        Store adaptive datum to file
        
        *Arguments*
            filename (string)
                Name of file to store in


    .. cpp:function:: void store(Table& table) const
    
        Store adaptive datum to :cpp:class:`Table`.
        
        *Arguments*
            table (:cpp:class:`Table`)
                Table to store in


    .. cpp:function:: void set_reference_value(const double reference)
    
        Set reference value for goal functional
        
        *Arguments*
            reference (double)
                The value.


