
.. Documentation for the header file dolfin/common/timing.h

.. _programmers_reference_cpp_common_timing:

timing.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void tic()

    Timing functions measure CPU time as determined by clock(),
    the precision of which seems to be 0.01 seconds.
    Start timing (should not be used internally in DOLFIN!)


.. cpp:function:: double toc()

    Return elapsed CPU time (should not be used internally in DOLFIN!)


.. cpp:function:: double time()

    Return current CPU time used by process


.. cpp:function:: Table timings(bool reset=false)

    Return a summary of timings and tasks in a Table, optionally clearing
    stored timings


.. cpp:function:: void list_timings(bool reset=false)

    List a summary of timings and tasks, optionally clearing stored timings


.. cpp:function:: void summary(bool reset=false)

    This function is deprecated, use list_timings


.. cpp:function:: double timing(std::string task, bool reset=false)

    Return timing (average) for given task, optionally clearing timing
    for task


