.. Documentation for the header file dolfin/mesh/MeshTopology.h

.. _programmers_reference_cpp_mesh_meshtopology:

MeshTopology.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MeshTopology

    MeshTopology stores the topology of a mesh, consisting of mesh entities
    and connectivity (incidence relations for the mesh entities). Note that
    the mesh entities don't need to be stored, only the number of entities
    and the connectivity. Any numbering scheme for the mesh entities is
    stored separately in a MeshFunction over the entities.
    
    A mesh entity e may be identified globally as a pair e = (dim, i), where
    dim is the topological dimension and i is the index of the entity within
    that topological dimension.

    .. cpp:function:: MeshTopology()
    
        Create empty mesh topology

    .. cpp:function:: MeshTopology(const MeshTopology& topology)
    
        Copy constructor

    .. cpp:function:: const MeshTopology& operator= (const MeshTopology& topology)
    
        Assignment

    .. cpp:function:: uint dim() const
    
        Return topological dimension

    .. cpp:function:: uint size(uint dim) const
    
        Return number of entities for given dimension

    .. cpp:function:: void clear()
    
        Clear all data

    .. cpp:function:: void init(uint dim)
    
        Initialize topology of given maximum dimension

    .. cpp:function:: void init(uint dim, uint size)
    
        Set number of entities (size) for given topological dimension

    .. cpp:function:: dolfin::MeshConnectivity& operator() (uint d0, uint d1)
    
        Return connectivity for given pair of topological dimensions

    .. cpp:function:: const dolfin::MeshConnectivity& operator() (uint d0, uint d1) const
    
        Return connectivity for given pair of topological dimensions

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

