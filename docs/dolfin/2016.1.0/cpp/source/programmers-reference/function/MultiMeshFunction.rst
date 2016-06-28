
.. Documentation for the header file dolfin/function/MultiMeshFunction.h

.. _programmers_reference_cpp_function_multimeshfunction:

MultiMeshFunction.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshFunction

    This class represents a function on a cut and composite finite
    element function space (MultiMesh) defined on one or more possibly
    intersecting meshes.


    .. cpp:function:: explicit MultiMeshFunction(std::shared_ptr<const MultiMeshFunctionSpace> V)
    
        Create MultiMesh function on given MultiMesh function space
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The MultiMesh function space.


    .. cpp:function:: std::shared_ptr<const Function> part(std::size_t i) const
    
        Return function (part) number i
        
        *Returns*
            :cpp:class:`Function`
                Function (part) number i


    .. cpp:function:: std::shared_ptr<GenericVector> vector()
    
        Return vector of expansion coefficients (non-const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients.


    .. cpp:function:: std::shared_ptr<const GenericVector> vector() const
    
        Return vector of expansion coefficients (const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients (const).


