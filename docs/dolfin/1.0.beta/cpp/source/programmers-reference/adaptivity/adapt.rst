
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


.. cpp:function:: const FunctionSpace& adapt(const FunctionSpace& space, boost::shared_ptr<const Mesh> refined_mesh)

    Refine function space based on refined mesh


.. cpp:function:: const Function& adapt(const Function& function, boost::shared_ptr<const Mesh> refined_mesh)

    Refine coefficient based on refined mesh


.. cpp:function:: const MeshFunction<dolfin::uint>& adapt(const MeshFunction<uint>& mesh_function, boost::shared_ptr<const Mesh> refined_mesh)

    Refine mesh function<uint> based on mesh


.. cpp:function:: const DirichletBC& adapt(const DirichletBC& bc, boost::shared_ptr<const Mesh> refined_mesh, const FunctionSpace& S)

    Refine Dirichlet bc based on refined mesh


.. cpp:function:: void adapt_markers(std::vector<std::pair<uint, uint> >& refined_markers, const Mesh& refined_mesh, const std::vector<std::pair<uint, uint> >& markers, const Mesh& mesh)

    Helper function for refinement of boundary conditions


.. cpp:function:: const Form& adapt(const Form& form, boost::shared_ptr<const Mesh> refined_mesh)

    Refine form based on refined mesh


.. cpp:function:: const LinearVariationalProblem& adapt(const LinearVariationalProblem& problem, boost::shared_ptr<const Mesh> refined_mesh)

    Refine linear variational problem based on mesh


.. cpp:function:: const NonlinearVariationalProblem& adapt(const NonlinearVariationalProblem& problem, boost::shared_ptr<const Mesh> refined_mesh)

    Refine nonlinear variational problem based on mesh


.. cpp:function:: const ErrorControl& adapt(const ErrorControl& ec, boost::shared_ptr<const Mesh> refined_mesh)

    Refine error control object based on mesh


