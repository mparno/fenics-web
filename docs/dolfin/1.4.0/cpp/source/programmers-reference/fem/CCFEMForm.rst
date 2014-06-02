
.. Documentation for the header file dolfin/fem/CCFEMForm.h

.. _programmers_reference_cpp_fem_ccfemform:

CCFEMForm.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMForm

    This class represents a variational form on a cut and composite
    finite element function space (CCFEM) defined on one or more
    possibly intersecting meshes.
    
    FIXME: Document usage of class with add() followed by build()


    .. cpp:function:: CCFEMForm(std::shared_ptr<const CCFEMFunctionSpace> function_space)
    
        Create empty linear CCFEM variational form (shared pointer version)


    .. cpp:function:: CCFEMForm(const CCFEMFunctionSpace& function_space)
    
        Create empty linear CCFEM variational form (reference version)


    .. cpp:function:: CCFEMForm(std::shared_ptr<const CCFEMFunctionSpace> function_space_0, std::shared_ptr<const CCFEMFunctionSpace> function_space_1)
    
        Create empty bilinear CCFEM variational form (shared pointer version)


    .. cpp:function:: CCFEMForm(const CCFEMFunctionSpace& function_space_0, const CCFEMFunctionSpace& function_space_1)
    
        Create empty bilinear CCFEM variational form (reference version)


    .. cpp:function:: std::size_t rank() const
    
        Return rank of form (bilinear form = 2, linear form = 1,
        functional = 0, etc)
        
        *Returns*
            std::size_t
                The rank of the form.


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number of forms (parts) of the CCFEM form
        
        *Returns*
            std::size_t
                The number of forms (parts) of the CCFEM form.


    .. cpp:function:: std::shared_ptr<const Form> part(std::size_t i) const
    
        Return form (part) number i
        
        *Returns*
            :cpp:class:`Form`
                Form (part) number i.


    .. cpp:function:: std::shared_ptr<const CCFEMFunctionSpace> function_space(std::size_t i) const
    
        Return function space for given argument
        
        *Arguments*
            i (std::size_t)
                Index
        
        *Returns*
            :cpp:class:`CCFEMFunctionSpace`
                Function space shared pointer.


    .. cpp:function:: void add(std::shared_ptr<const Form> form)
    
        Add form (shared pointer version)
        
        *Arguments*
            form (:cpp:class:`Form`)
                The form.


    .. cpp:function:: void add(const Form& form)
    
        Add form (reference version)
        
        *Arguments*
            form (:cpp:class:`Form`)
                The form.


    .. cpp:function:: void build()
    
        Build CCFEM form


    .. cpp:function:: void clear()
    
        Clear CCFEM form


