.. Documentation for the header file dolfin/log/log.h

.. _programmers_reference_cpp_log_log:

log.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: void info(std::string msg, ...)
    
        The DOLFIN log system provides the following set of functions for
        uniform handling of log messages, warnings and errors. In addition,
        macros are provided for debug messages and assertions.
        
        Only messages with a debug level higher than or equal to the current
        log level are printed (the default being zero). Logging may also be
        turned off by calling set_log_active(false).
        Print message

    .. cpp:function:: void info(const Parameters& parameters, bool verbose=false)
    
        Print parameter (using output of str() method)

    .. cpp:function:: void info(const Variable& variable, bool verbose=false)
    
        Print variable (using output of str() method)

    .. cpp:function:: void info_stream(std::ostream& out, std::string msg)
    
        Print message to stream

    .. cpp:function:: void info_underline(std:: string msg, ...)
    
        Print underlined message

    .. cpp:function:: void warning(std::string msg, ...)
    
        Print warning

    .. cpp:function:: void error(std::string msg, ...)
    
        Print error message and throw an exception

    .. cpp:function:: void log(int debug_level, std::string msg, ...)
    
        Print message at given debug level

    .. cpp:function:: void begin(std::string msg, ...)
    
        Begin task (increase indentation level)

    .. cpp:function:: void begin(int debug_level, std::string msg, ...)
    
        Begin task (increase indentation level)

    .. cpp:function:: void end()
    
        End task (decrease indentation level)

    .. cpp:function:: void set_log_active(bool active=true)
    
        Turn logging on or off (deprecated)

    .. cpp:function:: void logging(bool active=true)
    
        Turn logging on or off (deprecated, will be removed)

    .. cpp:function:: void set_log_level(int level)
    
        Set log level

    .. cpp:function:: void set_output_stream(std::ostream& out)
    
        Set output stream

    .. cpp:function:: int get_log_level()
    
        Get log level

    .. cpp:function:: void summary(bool reset=false)
    
        Print summary of timings and tasks, optionally clearing stored timings

    .. cpp:function:: double timing(std::string task, bool reset=false)
    
        Return timing (average) for given task, optionally clearing timing for task

    .. cpp:function:: void not_working_in_parallel(std::string what)
    
        Report that functionality has not (yet) been implemented to work in parallel

    .. cpp:function:: void check_equal(uint value, uint valid_value, std::string task, std::string value_name)
    
        Check value and print an informative error message if invalid

