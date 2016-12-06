
.. Documentation for the header file dolfin/mesh/MeshRenumbering.h

.. _programmers_reference_cpp_mesh_meshrenumbering:

MeshRenumbering.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshRenumbering

    This class implements renumbering algorithms for meshes.


    .. cpp:function:: static Mesh renumber_by_color(const Mesh& mesh, std::vector<std::size_t> coloring)
    
        Renumber mesh entities by coloring. This function is currently
        restricted to renumbering by cell coloring. The cells
        (cell-vertex connectivity) and the coordinates of the mesh are
        renumbered to improve the locality within each color. It is
        assumed that the mesh has already been colored and that only
        cell-vertex connectivity exists as part of the mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Mesh to be renumbered.
            coloring (_std::vector<std::size_t>_)
                Mesh coloring type.
        *Returns*
            :cpp:class:`Mesh`


