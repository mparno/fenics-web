
.. Documentation for the header file dolfin/log/log.h

.. _programmers_reference_cpp_log_log:

log.h
=====

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void info(std::string msg, ...)

    The DOLFIN log system provides the following set of functions for
    uniform handling of log messages, warnings and errors. In addition,
    macros are provided for debug messages and dolfin_assertions.
    
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


.. cpp:function:: void info_underline(std::string msg, ...)

    Print underlined message


.. cpp:function:: void warning(std::string msg, ...)

    Print warning


.. cpp:function:: void error(std::string msg, ...)

    Print error message and throw an exception.
    Note to developers: this function should not be used internally
    in DOLFIN. Use the more informative dolfin_error instead.


.. cpp:function:: void dolfin_error(std::string location, std::string task, std::string reason, ...)

    Print error message. Prefer this to the above generic error message.
    
    *Arguments*
        location (std::string)
            Name of the file from which the error message was generated.
        task (std::string)
            Name of the task that failed.
            Note that this string should begin with lowercase.
            Note that this string should not be punctuated.
        reason (std::string)
            A format string explaining the reason for the failure.
            Note that this string should begin with uppercase.
            Note that this string should not be punctuated.
            Note that this string may contain printf style formatting.
        ... (primitive types like int, std::size_t, double, bool)
            Optional arguments for the format string.
    
    Developers should read the file dolfin/log/README in the DOLFIN
    source tree for further notes about the use of this function.


.. cpp:function:: void deprecation(std::string feature, std::string version_deprecated, std::string message, ...)

    Issue deprecation warning for removed feature
    
    *Arguments*
        feature (std::string)
           Name of the feature that has been removed.
        version_deprecated (std::string)
           Version number of the release in which the feature is deprecated.
        message (std::string)
           A format string explaining the deprecation.


.. cpp:function:: void log(int debug_level, std::string msg, ...)

    Print message at given debug level


.. cpp:function:: void begin(std::string msg, ...)

    Begin task (increase indentation level)


.. cpp:function:: void begin(int debug_level, std::string msg, ...)

    Begin task (increase indentation level)


.. cpp:function:: void end()

    End task (decrease indentation level)


.. cpp:function:: void set_log_active(bool active=true)

    Turn logging on or off


.. cpp:function:: void set_log_level(int level)

    Set log level


.. cpp:function:: void set_output_stream(std::ostream& out)

    Set output stream


.. cpp:function:: int get_log_level()

    Get log level


.. cpp:function:: void monitor_memory_usage()

    Monitor memory usage. Call this function at the start of a
    program to continuously monitor the memory usage of the process.


.. cpp:function:: void not_working_in_parallel(std::string what)

    Report that functionality has not (yet) been implemented to work
    in parallel


