.. Documentation for the header file dolfin/common/Variable.h

.. _programmers_reference_cpp_common_variable:

Variable.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Variable

        Common base class for DOLFIN variables.

    .. cpp:function:: Variable()
    
        Create unnamed variable

    .. cpp:function:: Variable(const std::string name, const std::string label)
    
        Create variable with given name and label

    .. cpp:function:: const std::string& label() const
    
        Return label (description)

    .. cpp:function:: const std::string& name()  const
    
        Return name

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual ~Variable()
    
        Destructor

    .. cpp:function:: void disp() const
    
        Deprecated, to be removed

    .. cpp:function:: void rename(const std::string name, const std::string label)
    
        Rename variable

