
.. Documentation for the header file dolfin/generation/ImplicitDomainMeshGenerator.h

.. _programmers_reference_cpp_generation_implicitdomainmeshgenerator:

ImplicitDomainMeshGenerator.h
=============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ImplicitDomainMeshGenerator

    Generate surface or volume meshes from implicit representations
    of a domain modelled by the abstract dolfin::ImplicitSurface
    class. CGAL is used for the mesh generation.


    .. cpp:function:: static void generate(Mesh& mesh, const ImplicitSurface& surface, double cell_size)
    
        Create volume mesh from implicit surface representation


    .. cpp:function:: static void generate_surface(Mesh& mesh, const ImplicitSurface& surface, double cell_size)
    
        Create surface mesh from implicit surface representation. This
        function uses the CGAL 3D mesh genrator


