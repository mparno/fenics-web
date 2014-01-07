
.. Documentation for the header file dolfin/common/Hierarchical.h

.. _programmers_reference_cpp_common_hierarchical:

Hierarchical.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Hierarchical

    This class provides storage and data access for hierarchical
    classes; that is, classes where an object may have a child
    and a parent.
    
    Note to developers: each subclass of Hierarchical that
    implements an assignment operator must call the base class
    assignment operator at the *end* of the subclass assignment
    operator. See the Mesh class for an example.


    .. cpp:function:: Hierarchical(T& self)
    
        Constructor


    .. cpp:function:: std::size_t depth() const
    
        Return depth of the hierarchy; that is, the total number of
        objects in the hierarchy linked to the current object via
        child-parent relationships, including the object itself.
        
        *Returns*
            std::size_t
                The depth of the hierarchy.


    .. cpp:function:: bool has_parent() const
    
        Check if the object has a parent.
        
        *Returns*
            bool
                The return value is true iff the object has a parent.


    .. cpp:function:: bool has_child() const
    
        Check if the object has a child.
        
        *Returns*
            bool
                The return value is true iff the object has a child.


    .. cpp:function:: T& parent()
    
        Return parent in hierarchy. An error is thrown if the object
        has no parent.
        
        *Returns*
            _Object_
                The parent object.


    .. cpp:function:: const T& parent() const
    
        Return parent in hierarchy (const version).


    .. cpp:function:: boost::shared_ptr<T> parent_shared_ptr()
    
        Return shared pointer to parent. A zero pointer is returned if
        the object has no parent.
        
        *Returns*
            shared_ptr<T>
                The parent object.


    .. cpp:function:: boost::shared_ptr<const T> parent_shared_ptr() const
    
        Return shared pointer to parent (const version).


    .. cpp:function:: T& child()
    
        Return child in hierarchy. An error is thrown if the object
        has no child.
        
        *Returns*
            _T_
                The child object.


    .. cpp:function:: const T& child() const
    
        Return child in hierarchy (const version).


    .. cpp:function:: boost::shared_ptr<T> child_shared_ptr()
    
        Return shared pointer to child. A zero pointer is returned if
        the object has no child.
        
        *Returns*
            shared_ptr<T>
                The child object.


    .. cpp:function:: boost::shared_ptr<const T> child_shared_ptr() const
    
        Return shared pointer to child (const version).


    .. cpp:function:: T& root_node()
    
        Return root node object in hierarchy.
        
        *Returns*
            _T_
                The root node object.


    .. cpp:function:: const T& root_node() const
    
        Return root node object in hierarchy (const version).


    .. cpp:function:: boost::shared_ptr<T> root_node_shared_ptr()
    
        Return shared pointer to root node object in hierarchy.
        
        *Returns*
            _T_
                The root node object.


    .. cpp:function:: boost::shared_ptr<const T> root_node_shared_ptr() const
    
        Return shared pointer to root node object in hierarchy (const version).


    .. cpp:function:: T& leaf_node()
    
        Return leaf node object in hierarchy.
        
        *Returns*
            _T_
                The leaf node object.


    .. cpp:function:: const T& leaf_node() const
    
        Return leaf node object in hierarchy (const version).


    .. cpp:function:: boost::shared_ptr<T> leaf_node_shared_ptr()
    
        Return shared pointer to leaf node object in hierarchy.
        
        *Returns*
            _T_
                The leaf node object.


    .. cpp:function:: boost::shared_ptr<const T> leaf_node_shared_ptr() const
    
        Return shared pointer to leaf node object in hierarchy (const version).


    .. cpp:function:: void set_parent(boost::shared_ptr<T> parent)
    
        Set parent


    .. cpp:function:: void clear_child()
    
        Clear child


    .. cpp:function:: void set_child(boost::shared_ptr<T> child)
    
        Set child


    .. cpp:function:: const Hierarchical& operator= (const Hierarchical& hierarchical)
    
        Assignment operator


    .. cpp:function:: void _debug() const
    
        Function useful for debugging the hierarchy


