
.. Documentation for the header file dolfin/mesh/MeshHierarchy.h

.. _programmers_reference_cpp_mesh_meshhierarchy:

MeshHierarchy.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshHierarchy

    .. cpp:function:: MeshHierarchy()
    
        Constructor


    .. cpp:function:: explicit MeshHierarchy(std::shared_ptr<const Mesh> mesh)
    
        Constructor with initial mesh


    .. cpp:function:: unsigned int size() const
    
        Number of meshes


    .. cpp:function:: std::shared_ptr<const Mesh> operator[](int i) const
    
        Get Mesh i, in range [0:size()] where 0 is the coarsest Mesh.


    .. cpp:function:: std::shared_ptr<const Mesh> finest() const
    
        Get the finest mesh of the MeshHierarchy


    .. cpp:function:: std::shared_ptr<const Mesh> coarsest() const
    
        Get the coarsest mesh of the MeshHierarchy


    .. cpp:function:: std::shared_ptr<const MeshHierarchy> refine (const MeshFunction<bool>& markers) const
    
        Refine finest mesh of existing hierarchy, creating a new hierarchy
        (level n -> n+1)


    .. cpp:function:: std::shared_ptr<const MeshHierarchy> unrefine() const
    
        Unrefine by returning the previous MeshHierarchy
        (level n -> n-1)
        Returns NULL for a MeshHierarchy containing a single Mesh


    .. cpp:function:: std::shared_ptr<const MeshHierarchy> coarsen (const MeshFunction<bool>& markers) const
    
        Coarsen finest mesh by one level, based on markers (level n->n)


    .. cpp:function:: std::vector<std::size_t> weight() const
    
        Calculate the number of cells on the finest Mesh
        which are descendents of each cell on the coarsest Mesh,
        returning a vector over the cells of the coarsest Mesh.


    .. cpp:function:: std::shared_ptr<Mesh> rebalance() const
    
        Rebalance across processes


