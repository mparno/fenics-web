
.. Documentation for the header file dolfin/fem/MultiMeshForm.h

.. _programmers_reference_cpp_fem_multimeshform:

MultiMeshForm.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshForm

    This class represents a variational form on a cut and composite
    finite element function space (MultiMesh) defined on one or more
    possibly intersecting meshes.


    .. cpp:function:: MultiMeshForm(std::shared_ptr<const MultiMesh> multimesh)
    
        Create empty multimesh functional


    .. cpp:function:: MultiMeshForm(std::shared_ptr<const MultiMeshFunctionSpace> function_space)
    
        Create empty linear multimesh variational form


    .. cpp:function:: MultiMeshForm(std::shared_ptr<const MultiMeshFunctionSpace> function_space_0, std::shared_ptr<const MultiMeshFunctionSpace> function_space_1)
    
        Create empty bilinear multimesh variational form


    .. cpp:function:: std::size_t rank() const
    
        Return rank of form (bilinear form = 2, linear form = 1,
        functional = 0, etc)
        
        *Returns*
            std::size_t
                The rank of the form.


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number of forms (parts) of the MultiMesh form
        
        *Returns*
            std::size_t
                The number of forms (parts) of the MultiMesh form.


    .. cpp:function:: std::shared_ptr<const MultiMesh> multimesh() const
    
        Extract common multimesh from form
        
        *Returns*
            :cpp:class:`MultiMesh`
                The mesh.


    .. cpp:function:: std::shared_ptr<const Form> part(std::size_t i) const
    
        Return form (part) number i
        
        *Returns*
            :cpp:class:`Form`
                Form (part) number i.


    .. cpp:function:: std::shared_ptr<const MultiMeshFunctionSpace> function_space(std::size_t i) const
    
        Return function space for given argument
        
        *Arguments*
            i (std::size_t)
                Index
        
        *Returns*
            :cpp:class:`MultiMeshFunctionSpace`
                Function space shared pointer.


    .. cpp:function:: void add(std::shared_ptr<const Form> form)
    
        Add form (shared pointer version)
        
        *Arguments*
            form (:cpp:class:`Form`)
                The form.


    .. cpp:function:: void build()
    
        Build MultiMesh form


    .. cpp:function:: void clear()
    
        Clear MultiMesh form


