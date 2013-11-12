
.. Documentation for the header file dolfin/generation/EllipseMesh.h

.. _programmers_reference_cpp_generation_ellipsemesh:

EllipseMesh.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EllipseMesh

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Triangular mesh of an ellipse. CGAL is used to generate the
    mesh.


    .. cpp:function:: EllipseMesh(Point p, std::vector<double> ellipse_dims, double cell_size)
    
                EllipsoidMesh mesh(Point(1.0, 2.0), dims, 0.2);
        


.. cpp:class:: CircleMesh

    *Parent class(es)*
    
        * :cpp:class:`EllipseMesh`
        
    .. cpp:function:: CircleMesh(Point p, double radius, double cell_size)
    
        Create an unstructured :cpp:class:`Mesh` of a circle
        
        *Arguments*
            center (:cpp:class:`Point`)
                Center of the ellipsoid
            radius (double)
                Axes lengths
            cell_size (double)
                Cell size measure
        
        *Example*
            .. code-block:: c++
        
                // Create sphere with center (1.0, 2.0) and
                // radius 3.0. Cell size is 0.2.
                SphereMesh mesh(Point(1.0, 2.0), 3.0, 0.2);
        


