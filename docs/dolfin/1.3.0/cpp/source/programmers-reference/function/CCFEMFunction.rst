
.. Documentation for the header file dolfin/function/CCFEMFunction.h

.. _programmers_reference_cpp_function_ccfemfunction:

CCFEMFunction.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMFunction

    This class represents a function on a cut and composite finite
    element function space (CCFEM) defined on one or more possibly
    intersecting meshes.


    .. cpp:function:: explicit CCFEMFunction(const CCFEMFunctionSpace& V)
    
        Create CCFEM function on given CCFEM function space
        
        *Arguments*
            V (:cpp:class:`CCFEMFunctionSpace`)
                The CCFEM function space.
        
        *Example*
            .. code-block:: c++
        
                CCFEMFunction u(V);
        


    .. cpp:function:: explicit CCFEMFunction(boost::shared_ptr<const CCFEMFunctionSpace> V)
    
        Create CCFEM function on given CCFEM function space (shared
        pointer version)
        
        *Arguments*
            V (:cpp:class:`CCFEMFunctionSpace`)
                The CCFEM function space.


    .. cpp:function:: const Function& part(std::size_t i) const
    
        Return function (part) number i
        
        *Returns*
            :cpp:class:`Function`
                Function (part) number i


    .. cpp:function:: boost::shared_ptr<GenericVector> vector()
    
        Return vector of expansion coefficients (non-const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients.


    .. cpp:function:: boost::shared_ptr<const GenericVector> vector() const
    
        Return vector of expansion coefficients (const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients (const).


