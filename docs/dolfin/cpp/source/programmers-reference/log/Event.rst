
.. Documentation for the header file dolfin/log/Event.h

.. _programmers_reference_cpp_log_event:

Event.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Event

    A event is a string message which is displayed
    only a limited number of times.
    
    *Example*
        .. code-block:: c++
    
            Event event("System is stiff, damping is needed.");
            while ()
            {
              ...
              if ( ... )
              {
                event();
                ...
              }
            }


    .. cpp:function:: Event(const std::string msg, unsigned int maxcount = 1)
    
        Constructor


    .. cpp:function:: void operator() ()
    
        Display message


    .. cpp:function:: unsigned int count() const
    
        Display count


    .. cpp:function:: unsigned int maxcount() const
    
        Maximum display count


