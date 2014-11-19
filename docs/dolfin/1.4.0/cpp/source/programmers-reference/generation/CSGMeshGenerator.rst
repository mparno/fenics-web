
.. Documentation for the header file dolfin/generation/CSGMeshGenerator.h

.. _programmers_reference_cpp_generation_csgmeshgenerator:

CSGMeshGenerator.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CSGMeshGenerator

    Mesh generator for Constructive Solid Geometry (CSG)


    .. cpp:function:: static void generate(Mesh& mesh, const CSGGeometry& geometry, std::size_t resolution)
    
        Generate mesh from CSG geometry


    .. cpp:function:: static void generate(BoundaryMesh& mesh, const CSGGeometry& geometry)
    
        Generate boundary mesh from the surface of a CSG geometry


