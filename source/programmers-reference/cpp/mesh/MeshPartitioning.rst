.. Documentation for the header file dolfin/mesh/MeshPartitioning.h

.. _programmers_reference_cpp_mesh_meshpartitioning:

MeshPartitioning.h
==================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: MeshPartitioning

    This class partitions and distributes a mesh based on
    partitioned local mesh data. Note that the local mesh data will
    also be repartitioned and redistributed during the computation
    of the mesh partitioning.
    
    After partitioning, each process has a local mesh and set of
    mesh data that couples the meshes together.
    
    The following mesh data is created:
    
    1. "global entity indices 0" (MeshFunction<uint>)
    
    This maps each local vertex to its global index.
    
    2. "overlap" (std::map<uint, std::vector<uint> >)
    
    This maps each shared vertex to a list of the processes sharing
    the vertex.
    
    3. "global entity indices %d" (MeshFunction<uint>)
    
    After partitioning, the function number_entities() may be called
    to create global indices for all entities of a given topological
    dimension. These are stored as mesh data (MeshFunction<uint>)
    named
    
       "global entity indices 1"
       "global entity indices 2"
       etc
    
    4. "num global entities" (std::vector<uint>)
    
    The function number_entities also records the number of global
    entities for the dimension of the numbered entities in the array
    named "num global entities". This array has size D + 1, where D
    is the topological dimension of the mesh. This array is
    initially created by the mesh and then contains only the number
    entities of dimension 0 (vertices) and dimension D (cells).

    .. cpp:function:: static void number_entities(const Mesh& mesh, uint d)
    
        Create global entity indices for entities of dimension d

    .. cpp:function:: static void partition(Mesh& mesh)
    
        Create a partitioned mesh based on local meshes

    .. cpp:function:: static void partition(Mesh& mesh, LocalMeshData& data)
    
        Create a partitioned mesh based on local mesh data

