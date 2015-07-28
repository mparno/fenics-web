
.. Documentation for the header file dolfin/adaptivity/adapt.h

.. _programmers_reference_cpp_adaptivity_adapt:

adapt.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: const Mesh& adapt(const Mesh& mesh)

    Refine mesh uniformly


.. cpp:function:: const Mesh& adapt(const Mesh& mesh, const MeshFunction<bool>& cell_markers)

    Refine mesh based on cell markers


.. cpp:function:: const FunctionSpace& adapt(const FunctionSpace& space)

    Refine function space uniformly


.. cpp:function:: const FunctionSpace& adapt(const FunctionSpace& space, const MeshFunction<bool>& cell_markers)

    Refine function space based on cell markers


.. cpp:function:: const FunctionSpace& adapt(const FunctionSpace& space, std::shared_ptr<const Mesh> adapted_mesh)

    Refine function space based on refined mesh


.. cpp:function:: const Function& adapt(const Function& function, std::shared_ptr<const Mesh> adapted_mesh, bool interpolate=true)

    Adapt Function based on adapted mesh
    
    *Arguments*
        function  (:cpp:class:`Function`)
            The function that should be adapted
        adapted_mesh  (:cpp:class:`Mesh`)
            The new mesh
        interpolate (bool)
            Optional argument, default is true. If false, the
            function's function space is adapted, but the values are
            not interpolated.
    
    *Returns*
        _Function__
            The adapted function


.. cpp:function:: const GenericFunction& adapt(const GenericFunction& function, std::shared_ptr<const Mesh> adapted_mesh)

    Refine GenericFunction based on refined mesh


.. cpp:function:: const MeshFunction<std::size_t>& adapt(const MeshFunction<std::size_t>& mesh_function, std::shared_ptr<const Mesh> adapted_mesh)

    Refine mesh function<std::size_t> based on mesh


.. cpp:function:: const DirichletBC& adapt(const DirichletBC& bc, std::shared_ptr<const Mesh> adapted_mesh, const FunctionSpace& S)

    Refine Dirichlet bc based on refined mesh


.. cpp:function:: void adapt_markers(std::vector<std::size_t>& refined_markers, const Mesh& adapted_mesh, const std::vector<std::size_t>& markers, const Mesh& mesh)

    Helper function for refinement of boundary conditions


.. cpp:function:: const Form& adapt(const Form& form, std::shared_ptr<const Mesh> adapted_mesh, bool adapt_coefficients=true)

    Adapt form based on adapted mesh
    
    *Arguments*
        form  (:cpp:class:`Form`)
            The form that should be adapted
        adapted_mesh  (:cpp:class:`Mesh`)
            The new mesh
        adapt_coefficients (bool)
            Optional argument, default is true. If false, the form
            coefficients are not explicitly adapted, but pre-adapted
            coefficients will be transferred.
    
    *Returns*
        _Form__
            The adapted form


.. cpp:function:: const LinearVariationalProblem& adapt(const LinearVariationalProblem& problem, std::shared_ptr<const Mesh> adapted_mesh)

    Refine linear variational problem based on mesh


.. cpp:function:: const NonlinearVariationalProblem& adapt(const NonlinearVariationalProblem& problem, std::shared_ptr<const Mesh> adapted_mesh)

    Refine nonlinear variational problem based on mesh


.. cpp:function:: const ErrorControl& adapt(const ErrorControl& ec, std::shared_ptr<const Mesh> adapted_mesh, bool adapt_coefficients=true)

    Adapt error control object based on adapted mesh
    
    *Arguments*
        ec  (:cpp:class:`ErrorControl`)
            The error control object to be adapted
        adapted_mesh  (:cpp:class:`Mesh`)
            The new mesh
        adapt_coefficients (bool)
            Optional argument, default is true. If false, any form
            coefficients are not explicitly adapted, but pre-adapted
            coefficients will be transferred.
    
    *Returns*
        _ErrorControl__
            The adapted error control object


