.. Documentation for the header file dolfin/log/Progress.h

.. _programmers_reference_cpp_log_progress:

Progress.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Progress

    This class provides a simple way to create and update progress
    bars during a computation. A progress bar may be used either
    in an iteration with a known number of steps:
    
        Progress p("Iterating...", n);
        for (int i = 0; i < n; i++)
        {
          ...
          p++;
        }
    
    or in an iteration with an unknown number of steps:
    
        Progress p("Iterating...");
        while (t < T)
        {
          ...
          p = t / T;
        }

    .. cpp:function:: Progress(std::string title)
    
        Create progress bar with an unknown number of steps

    .. cpp:function:: Progress(std::string title, unsigned int n)
    
        Create progress bar with a known number of steps

    .. cpp:function:: void operator++(int)
    
        Increment progress

    .. cpp:function:: void operator=(double p)
    
        Set current position

    .. cpp:function:: ~Progress()
    
        Destructor

