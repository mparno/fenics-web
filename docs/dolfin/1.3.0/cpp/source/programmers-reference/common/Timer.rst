
.. Documentation for the header file dolfin/common/Timer.h

.. _programmers_reference_cpp_common_timer:

Timer.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Timer

    A timer can be used for timing tasks. The basic usage is
    
      Timer timer("Assembling over cells");
    
    The timer is started at construction and timing ends
    when the timer is destroyed (goes out of scope). It is
    also possible to start and stop a timer explicitly by
    
      timer.start();
      timer.stop();
    
    Timings are stored globally and a summary may be printed
    by calling
    
      list_timings();


    .. cpp:function:: Timer(std::string task)
    
        Create timer


    .. cpp:function:: void start()
    
        Start timer


    .. cpp:function:: double stop()
    
        Stop timer


    .. cpp:function:: double value() const
    
        Return value of timer (or time at start if not stopped)


