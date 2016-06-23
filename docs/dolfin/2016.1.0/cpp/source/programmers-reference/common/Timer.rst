
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


    .. cpp:function:: Timer()
    
        Create timer without logging


    .. cpp:function:: Timer(std::string task)
    
        Create timer with logging


    .. cpp:function:: void start()
    
        Zero and start timer


    .. cpp:function:: void resume()
    
        Resume timer. Not well-defined for logging timer


    .. cpp:function:: double stop()
    
        Stop timer, return wall time elapsed and store timing data
        into logger


    .. cpp:function:: std::tuple<double, double, double> elapsed() const
    
        Return wall, user and system time in seconds. Wall-clock time
        has precision around 1 microsecond; user and system around
        10 millisecond.


