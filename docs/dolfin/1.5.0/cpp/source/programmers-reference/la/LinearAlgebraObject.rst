
.. Documentation for the header file dolfin/la/LinearAlgebraObject.h

.. _programmers_reference_cpp_la_linearalgebraobject:

LinearAlgebraObject.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LinearAlgebraObject

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This is a common base class for all DOLFIN linear algebra
    objects. In particular, it provides casting mechanisms between
    different types.


    .. cpp:function:: const T& down_cast() const
    
        Cast object to its derived class, if possible (const version).
        An error is thrown if the cast is unsuccessful.


    .. cpp:function:: T& down_cast()
    
        Cast object to its derived class, if possible (non-const version).
        An error is thrown if the cast is unsuccessful.


    .. cpp:function:: static std::shared_ptr<X> down_cast(std::shared_ptr<Y> A)
    
        Cast shared pointer object to its derived class, if possible.
        Caller must check for success (returns null if cast fails).


    .. cpp:function:: const LinearAlgebraObject* instance() const
    
        Return concrete instance / unwrap (const version)


    .. cpp:function:: LinearAlgebraObject* instance()
    
        Return concrete instance / unwrap (non-const version)


    .. cpp:function:: std::shared_ptr<const LinearAlgebraObject> shared_instance() const
    
        Return concrete shared ptr instance / unwrap (const version)


    .. cpp:function:: std::shared_ptr<LinearAlgebraObject> shared_instance()
    
        Return concrete shared ptr instance / unwrap (non-const version)


    .. cpp:function:: Y& as_type(X& x)
    
        Cast object to its derived class, if possible (non-const version).
        An error is thrown if the cast is unsuccessful.


    .. cpp:function:: std::shared_ptr<Y> as_type(std::shared_ptr<X> x)
    
        Cast shared pointer object to its derived class, if possible.
        Caller must check for success (returns null if cast fails).


    .. cpp:function:: bool has_type(const X& x)
    
        Check whether the object matches a specific type


