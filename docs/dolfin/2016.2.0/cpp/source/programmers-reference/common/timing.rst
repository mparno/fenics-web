
.. Documentation for the header file dolfin/common/timing.h

.. _programmers_reference_cpp_common_timing:

timing.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: enum class TimingClear

    Parameter specifying whether to clear timing(s):
      * ``TimingClear::keep``
      * ``TimingClear::clear``


.. cpp:function:: enum class TimingType

    Timing types:
      * ``TimingType::wall`` wall-clock time
      * ``TimingType::user`` user (cpu) time
      * ``TimingType::system`` system (kernel) time
    
    Precision of wall is around 1 microsecond, user and system are around
    10 millisecond (on Linux).


.. cpp:function:: void tic()

    Start timing (should not be used internally in DOLFIN!)


.. cpp:function:: double toc()

    Return elapsed wall time (should not be used internally in DOLFIN!)


.. cpp:function:: double time()

    Return wall time elapsed since some implementation dependent epoch


.. cpp:function:: Table timings(TimingClear clear, std::set<TimingType> type)

    Return a summary of timings and tasks in a :cpp:class:`Table`, optionally clearing
    stored timings
    
    *Arguments*
        clear (TimingClear)
            * ``TimingClear::clear`` resets stored timings
            * ``TimingClear::keep`` leaves stored timings intact
        type (std::set<TimingType>)
            subset of ``{ TimingType::wall, TimingType::user,
            TimingType::system }``
    
    *Returns*
        :cpp:class:`Table`
            :cpp:class:`Table` with timings


.. cpp:function:: void list_timings(TimingClear clear, std::set<TimingType> type)

    List a summary of timings and tasks, optionally clearing stored timings.
    ``MPI_AVG`` reduction is printed. Collective on ``MPI_COMM_WORLD``.
    
    *Arguments*
        clear (TimingClear)
            * ``TimingClear::clear`` resets stored timings
            * ``TimingClear::keep`` leaves stored timings intact
        type (std::set<TimingType>)
            subset of ``{ TimingType::wall, TimingType::user,
            TimingType::system }``


.. cpp:function:: void dump_timings_to_xml(std::string filename, TimingClear clear)

    Dump a summary of timings and tasks to XML file, optionally clearing
    stored timings. ``MPI_MAX``, ``MPI_MIN`` and ``MPI_AVG`` reductions are
    stored. Collective on ``MPI_COMM_WORLD``.
    
    *Arguments*
        filename (std::string)
            output filename; must have ``.xml`` suffix; existing file
            is silently overwritten
        clear (TimingClear)
            * ``TimingClear::clear`` resets stored timings
            * ``TimingClear::keep`` leaves stored timings intact


.. cpp:function:: std::tuple<std::size_t, double, double, double> timing(std::string task, TimingClear clear)

    Return timing (count, total wall time, total user time,
    total system time) for given task, optionally clearing
    all timings for the task
    
    *Arguments*
        task (std::string)
            name of a task
        clear (TimingClear)
            * ``TimingClear::clear`` resets stored timings
            * ``TimingClear::keep`` leaves stored timings intact
    
    *Returns*
        std::tuple<std::size_t, double, double, double>
            (count, total wall time, total user time, total system time)


