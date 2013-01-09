
.. Documentation for the header file dolfin/mesh/MeshPartitioning.h

.. _programmers_reference_cpp_mesh_meshpartitioning:

MeshPartitioning.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshPartitioning

    This class partitions and distributes a mesh based on
    partitioned local mesh data.The local mesh data will
    also be repartitioned and redistributed during the computation
    of the mesh partitioning.
    
    After partitioning, each process has a local mesh and set of
    mesh data that couples the meshes together.
    
    3. "global entity indices %d" (MeshFunction<std::size_t>)
    
    After partitioning, the function number_entities() may be called
    to create global indices for all entities of a given topological
    dimension. These are stored as mesh data (MeshFunction<std::size_t>)
    named
    
       "global entity indices 1"
       "global entity indices 2"
       etc
    
    4. "num global entities" (std::vector<std::size_t>)
    
    The function number_entities also records the number of global
    entities for the dimension of the numbered entities in the array
    named "num global entities". This array has size D + 1, where D
    is the topological dimension of the mesh. This array is
    initially created by the mesh and then contains only the number
    entities of dimension 0 (vertices) and dimension D (cells).


    .. cpp:function:: static void build_distributed_mesh(Mesh& mesh)
    
        Build a partitioned mesh based on a local mesh


    .. cpp:function:: static void build_distributed_mesh(Mesh& mesh, const LocalMeshData& data)
    
        Build a partitioned mesh based on local mesh data


    .. cpp:function:: static void build_distributed_value_collection(MeshValueCollection<T>& values, const LocalMeshValueCollection<T>& local_data, const Mesh& mesh)
    
        Build a MeshValueCollection based on LocalMeshValueCollection


