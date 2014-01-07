
.. Documentation for the header file dolfin/generation/CSGOperators.h

.. _programmers_reference_cpp_generation_csgoperators:

CSGOperators.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CSGOperator

    *Parent class(es)*
    
        * :cpp:class:`CSGGeometry`
        
.. cpp:class:: CSGUnion

    *Parent class(es)*
    
        * :cpp:class:`CSGOperator`
        
    Union of CSG geometries


    .. cpp:function:: CSGUnion(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create union of two geometries


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: CSGDifference

    *Parent class(es)*
    
        * :cpp:class:`CSGOperator`
        
    Difference of CSG geometries


    .. cpp:function:: CSGDifference(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create union of two geometries


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: CSGIntersection

    *Parent class(es)*
    
        * :cpp:class:`CSGOperator`
        
    Intersection of CSG geometries


    .. cpp:function:: CSGIntersection(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create intersection of two geometries


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


    .. cpp:function:: boost::shared_ptr<CSGUnion> operator+(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create union of two geometries


    .. cpp:function:: boost::shared_ptr<CSGUnion> operator+(CSGGeometry& g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create union of two geometries


    .. cpp:function:: boost::shared_ptr<CSGUnion> operator+(boost::shared_ptr<CSGGeometry> g0, CSGGeometry& g1)
    
        Create union of two geometries


    .. cpp:function:: boost::shared_ptr<CSGUnion> operator+(CSGGeometry& g0, CSGGeometry& g1)
    
        Create union of two geometries


    .. cpp:function:: boost::shared_ptr<CSGDifference> operator-(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create difference of two geometries


    .. cpp:function:: boost::shared_ptr<CSGDifference> operator-(CSGGeometry& g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create difference of two geometries


    .. cpp:function:: boost::shared_ptr<CSGDifference> operator-(boost::shared_ptr<CSGGeometry> g0, CSGGeometry& g1)
    
        Create union of two geometries


    .. cpp:function:: boost::shared_ptr<CSGDifference> operator-(CSGGeometry& g0, CSGGeometry& g1)
    
        Create difference of two geometries


    .. cpp:function:: boost::shared_ptr<CSGIntersection> operator*(boost::shared_ptr<CSGGeometry> g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create intersection  of two geometries


    .. cpp:function:: boost::shared_ptr<CSGIntersection> operator*(CSGGeometry& g0, boost::shared_ptr<CSGGeometry> g1)
    
        Create intersection of two geometries


    .. cpp:function:: boost::shared_ptr<CSGIntersection> operator*(boost::shared_ptr<CSGGeometry> g0, CSGGeometry& g1)
    
        Create intersection of two geometries


    .. cpp:function:: boost::shared_ptr<CSGIntersection> operator*(CSGGeometry& g0, CSGGeometry& g1)
    
        Create intersection of two geometries


