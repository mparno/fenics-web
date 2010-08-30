.. Documentation for the header file dolfin/mesh/refine.h

.. _programmers_reference_cpp_mesh_refine:

refine.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: Mesh refine(const Mesh& mesh)
    
        Create uniformly refined mesh

    .. cpp:function:: Mesh refine(const Mesh& mesh, const MeshFunction<bool>& cell_markers)
    
        Create locally refined mesh

    .. cpp:function:: void refine(Mesh& refined_mesh, const Mesh& mesh)
    
        Create uniformly refined mesh

    .. cpp:function:: void refine(Mesh& refined_mesh, const Mesh& mesh, const MeshFunction<bool>& cell_markers)
    
        Create locally refined mesh

