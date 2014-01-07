
.. Documentation for the header file dolfin/generation/CSGGeometry.h

.. _programmers_reference_cpp_generation_csggeometry:

CSGGeometry.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CSGGeometry

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    Geometry described by Constructive Solid Geometry (CSG)


    .. cpp:function:: CSGGeometry()
    
        Constructor


    .. cpp:function:: std::size_t dim() const = 0
    
        Return dimension of geometry


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Informal string representation


    .. cpp:function:: void set_subdomain(std::size_t i, boost::shared_ptr<CSGGeometry> s)
    
        Define subdomain. This feature is 2D only.
        The subdomain is itself a CSGGeometry and the corresponding
        cells in the resulting will be marked with i
        If subdomains overlap, the latest added will take precedence.


