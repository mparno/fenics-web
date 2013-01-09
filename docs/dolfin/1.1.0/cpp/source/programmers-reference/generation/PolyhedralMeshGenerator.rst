
.. Documentation for the header file dolfin/generation/PolyhedralMeshGenerator.h

.. _programmers_reference_cpp_generation_polyhedralmeshgenerator:

PolyhedralMeshGenerator.h
=========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PolyhedralMeshGenerator

    Polyhedral mesh generator that uses CGAL


    .. cpp:function:: static void generate(Mesh& mesh, const std::string off_file, double cell_size, bool detect_sharp_features=true)
    
        Create mesh from Object File Format (.off) file


    .. cpp:function:: static void generate(Mesh& mesh, const std::vector<Point>& vertices, const std::vector<std::vector<std::size_t> >& facets, double cell_size, bool detect_sharp_features=true)
    
        Create mesh from a collection of facets


    .. cpp:function:: static void generate_surface_mesh(Mesh& mesh, const std::string off_file, double cell_size, bool detect_sharp_features=true)
    
        Create a surface mesh from Object File Format (.off) file


    .. cpp:function:: static void cgal_generate(Mesh& mesh, T& p, double cell_size, bool detect_sharp_features)
    
        Create mesh from a CGAL polyhedron


    .. cpp:function:: static void cgal_generate_surface_mesh(Mesh& mesh, T& p, double cell_size, bool detect_sharp_features)
    
        Create surface mesh from a CGAL polyhedron


