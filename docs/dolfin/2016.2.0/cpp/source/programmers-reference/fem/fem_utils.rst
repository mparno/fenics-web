
.. Documentation for the header file dolfin/fem/fem_utils.h

.. _programmers_reference_cpp_fem_fem_utils:

fem_utils.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: std::vector<std::size_t> dof_to_vertex_map(const FunctionSpace& space)

    Return a map between dofs indices and vertex indices
    
    Only works for FunctionSpace with dofs exclusively on vertices.
    For mixed FunctionSpaces vertex index is offset with the number
    of dofs per vertex.
    
    In parallel the returned map maps both owned and unowned dofs
    (using local indices) thus covering all the vertices. Hence the
    returned map is an inversion of _vertex_to_dof_map_.
    
    *Arguments*
        space (:cpp:class:`FunctionSpace`)
            The FunctionSpace for what the dof to vertex map should
            be computed for
    
    *Returns*
        std::vector<std::size_t>
            The dof to vertex map


.. cpp:function:: std::vector<dolfin::la_index> vertex_to_dof_map(const FunctionSpace& space)

    Return a map between vertex indices and dofs indices
    
    Only works for FunctionSpace with dofs exclusively on vertices.
    For mixed FunctionSpaces dof index is offset with the number of
    dofs per vertex.
    
    *Arguments*
        space (:cpp:class:`FunctionSpace`)
            The FunctionSpace for what the vertex to dof map should
            be computed for
    
    *Returns*
        std::vector<dolfin::la_index>
            The vertex to dof map


.. cpp:function:: void set_coordinates(MeshGeometry& geometry, const Function& position)

    Sets mesh coordinates from function
    
    Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
    (where d is topological dimension of the mesh and r is maximal
    dimension of entity associated with any coordinate node). Consider
    clearing unneeded connectivities when finished.
    
    *Arguments*
        geometry (:cpp:class:`MeshGeometry`)
            Mesh geometry to be set
        position (:cpp:class:`Function`)
            Vectorial Lagrange function with matching degree and mesh


.. cpp:function:: void get_coordinates(Function& position, const MeshGeometry& geometry)

    Stores mesh coordinates into function
    
    Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
    (where d is topological dimension of the mesh and r is maximal
    dimension of entity associated with any coordinate node). Consider
    clearing unneeded connectivities when finished.
    
    *Arguments*
        position (:cpp:class:`Function`)
            Vectorial Lagrange function with matching degree and mesh
        geometry (:cpp:class:`MeshGeometry`)
            Mesh geometry to be stored


.. cpp:function:: Mesh create_mesh(Function& position)

    Creates mesh from coordinate function
    
    Topology is given by underlying mesh of the function space and
    geometry is given by function values. Hence resulting mesh
    geometry has a degree of the function space degree. Geometry of
    function mesh is ignored.
    
    Mesh connectivities d-0, d-1, ..., d-r are built on function mesh
    (where d is topological dimension of the mesh and r is maximal
    dimension of entity associated with any coordinate node). Consider
    clearing unneeded connectivities when finished.
    
    *Arguments*
        position (:cpp:class:`Function`)
            Vectorial Lagrange function with of any degree
    
    *Returns*
        :cpp:class:`Mesh`
            The mesh


