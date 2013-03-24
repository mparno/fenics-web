
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
    
    After partitioning, each process has a local mesh and some data
    that couples the meshes together.


    .. cpp:function:: static void build_distributed_mesh(Mesh& mesh)
    
        Build a partitioned mesh based on a local mesh on process 0


    .. cpp:function:: static void build_distributed_mesh(Mesh& mesh, const LocalMeshData& data)
    
        Build a partitioned mesh based on local mesh data that is
        distributed across processes


    .. cpp:function:: static void build_distributed_value_collection(MeshValueCollection<T>& values, const LocalMeshValueCollection<T>& local_data, const Mesh& mesh)
    
        Build a MeshValueCollection based on LocalMeshValueCollection


