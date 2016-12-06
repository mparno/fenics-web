
.. Documentation for the header file dolfin/function/MultiMeshFunction.h

.. _programmers_reference_cpp_function_multimeshfunction:

MultiMeshFunction.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshFunction

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class represents a function on a cut and composite finite
    element function space (MultiMesh) defined on one or more possibly
    intersecting meshes.


    .. cpp:function:: MultiMeshFunction()
    
        Constructor


    .. cpp:function:: explicit MultiMeshFunction(std::shared_ptr<const MultiMeshFunctionSpace> V)
    
        Create MultiMesh function on given MultiMesh function space
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The MultiMesh function space.


    .. cpp:function:: MultiMeshFunction(std::shared_ptr<const MultiMeshFunctionSpace> V, std::shared_ptr<GenericVector> x)
    
        Create MultiMesh function on given MultiMesh function space with a given vector
        (shared data)
        
        *Warning: This constructor is intended for internal library use only*
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The multimesh function space.
            x (:cpp:class:`GenericVector`)
                The vector.


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


    .. cpp:function:: std::shared_ptr<const MultiMeshFunctionSpace> function_space() const
    
        Return shared pointer to multi mesh function space
        
        *Returns*
            :cpp:class:`MultiMeshFunctionSpace`
                Return the shared pointer.


