
.. Documentation for the header file dolfin/common/Variable.h

.. _programmers_reference_cpp_common_variable:

Variable.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Variable

    Common base class for DOLFIN variables.


    .. cpp:function:: Variable()
    
        Create unnamed variable


    .. cpp:function:: Variable(const std::string name, const std::string label)
    
        Create variable with given name and label


    .. cpp:function:: Variable(const Variable& variable)
    
        Copy constructor


    .. cpp:function:: void rename(const std::string name, const std::string label)
    
        Rename variable


    .. cpp:function:: std::string name()  const
    
        Return name


    .. cpp:function:: std::string label() const
    
        Return label (description)


    .. cpp:function:: std::size_t id() const
    
        Get unique identifier.
        
        *Returns*
            _std::size_t_
                The unique integer identifier associated with the object.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


