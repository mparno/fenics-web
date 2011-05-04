.. Documentation for the header file dolfin/log/Progress.h

.. _programmers_reference_cpp_log_progress:

Progress.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Progress

    .. cpp:function:: Progress(std::string title, unsigned int n)
    
        Create progress bar with a known number of steps

    .. cpp:function:: Progress(std::string title)
    
        Create progress bar with an unknown number of steps

    .. cpp:function:: void operator=(double p)
    
        Set current position

    .. cpp:function:: void operator++(int)
    
        Increment progress

