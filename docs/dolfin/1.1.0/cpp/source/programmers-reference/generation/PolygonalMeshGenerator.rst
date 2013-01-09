
.. Documentation for the header file dolfin/generation/PolygonalMeshGenerator.h

.. _programmers_reference_cpp_generation_polygonalmeshgenerator:

PolygonalMeshGenerator.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PolygonalMeshGenerator

    Polygonal mesh generator that uses CGAL


    .. cpp:function:: static void generate(Mesh& mesh, const std::vector<Point>& vertices, double cell_size)
    
        Generate mesh of a polygonal domain described by domain vertices


    .. cpp:function:: static void generate(Mesh& mesh, const T& polygon, double cell_size)
    
        Generate mesh of a domain described by a CGAL polygon


