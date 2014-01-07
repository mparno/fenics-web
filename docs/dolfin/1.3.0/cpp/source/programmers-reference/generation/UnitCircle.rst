
.. Documentation for the header file dolfin/generation/UnitCircle.h

.. _programmers_reference_cpp_generation_unitcircle:

UnitCircle.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitCircle

    *Parent class(es)*
    
        * :cpp:class:`UnitCircleMesh`
        
    .. cpp:function:: UnitCircle(std::size_t n, std::string diagonal="crossed", std::string transformation="rotsumn")
    
        Create a uniform finite element :cpp:class:`Mesh` over the unit circle.
        This class is deprecated. Use :cpp:class:`UnitCircleMesh`.
        
        *Arguments*
            n (std::size_t)
                Resolution of the mesh.
            diagonal (std::string)
                Optional argument: A std::string indicating
                the direction of the diagonals.
            transformation (std::string)
                Optional argument: A std::string indicating
                the type of transformation used.


