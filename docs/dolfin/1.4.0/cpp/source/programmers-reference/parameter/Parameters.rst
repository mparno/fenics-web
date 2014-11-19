
.. Documentation for the header file dolfin/parameter/Parameters.h

.. _programmers_reference_cpp_parameter_parameters:

Parameters.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Parameters

    This class stores a set of parameters. Each parameter is
    identified by a unique string (the key) and a value of some
    given value type. Parameter sets can be nested at arbitrary
    depths.
    
    A parameter may be either int, double, string or boolean valued.
    
    Parameters may be added as follows:
    
      Parameters p("my_parameters");
      p.add("relative_tolerance",  1e-15);
      p.add("absolute_tolerance",  1e-15);
      p.add("gmres_restart",       30);
      p.add("monitor_convergence", false);
    
    Parameters may be changed as follows:
    
      p["gmres_restart"] = 50;
    
    Parameter values may be retrieved as follows:
    
      int gmres_restart = p["gmres_restart"];
    
    Parameter sets may be nested as follows:
    
      Parameters q("nested_parameters");
      p.add(q);
    
    Nested parameters may then be accessed by
    
      p("nested_parameters")["..."]
    
    Parameters may be nested at arbitrary depths.
    
    Parameters may be parsed from the command-line as follows:
    
      p.parse(argc, argv);
    
    Note: spaces in parameter keys are not allowed (to simplify
    usage from command-line).


    .. cpp:function:: explicit Parameters(std::string key = "parameters")
    
        Create empty parameter set


    .. cpp:function:: Parameters(const Parameters& parameters)
    
        Copy constructor


    .. cpp:function:: std::string name() const
    
        Return name for parameter set


    .. cpp:function:: void rename(std::string key)
    
        Rename parameter set


    .. cpp:function:: void clear()
    
        Clear parameter set


    .. cpp:function:: void add(std::string key)
    
        Add an unset parameter of type T. For example, to create a
        unset parameter of type bool, do
        parameters.add<bool>("my_setting")


    .. cpp:function:: void add(std::string key, T min, T max)
    
        Add an unset parameter of type T with allows parameters. For
        example, to create a unset parameter of type bool, do
        parameters.add<bool>("my_setting")


    .. cpp:function:: void add(std::string key, std::set<T> valid_values)
    
        Add an unset parameter of type T with allows parameters. For
        example, to create a unset parameter of type bool, do
        parameters.add<bool>("my_setting")


    .. cpp:function:: void add(std::string key, int value)
    
        Add int-valued parameter


    .. cpp:function:: void add(std::string key, int value, int min_value, int max_value)
    
        Add int-valued parameter with given range


    .. cpp:function:: void add(std::string key, double value)
    
        Add double-valued parameter


    .. cpp:function:: void add(std::string key, double value, double min_value, double max_value)
    
        Add double-valued parameter with given range


    .. cpp:function:: void add(std::string key, std::string value)
    
        Add string-valued parameter


    .. cpp:function:: void add(std::string key, const char* value)
    
        Add string-valued parameter


    .. cpp:function:: void add(std::string key, std::string value, std::set<std::string> range)
    
        Add string-valued parameter with given range


    .. cpp:function:: void add(std::string key, const char* value, std::set<std::string> range)
    
        Add string-valued parameter with given range


    .. cpp:function:: void add(std::string key, bool value)
    
        Add bool-valued parameter


    .. cpp:function:: void add(const Parameters& parameters)
    
        Add nested parameter set


    .. cpp:function:: void remove(std::string key)
    
        Remove parameter or parameter set with given key


    .. cpp:function:: void parse(int argc, char* argv[])
    
        Parse parameters from command-line


    .. cpp:function:: void update(const Parameters& parameters)
    
        Update parameters with another set of parameters


    .. cpp:function:: Parameter& operator[] (std::string key)
    
        Return parameter for given key


    .. cpp:function:: const Parameter& operator[] (std::string key) const
    
        Return parameter for given key (const version)


    .. cpp:function:: Parameters& operator() (std::string key)
    
        Return nested parameter set for given key


    .. cpp:function:: const Parameters& operator() (std::string key) const
    
        Return nested parameter set for given key (const)


    .. cpp:function:: const Parameters& operator= (const Parameters& parameters)
    
        Assignment operator


    .. cpp:function:: bool has_key(std::string key) const
    
        Check if parameter set has key (parameter or nested parameter set)


    .. cpp:function:: bool has_parameter(std::string key) const
    
        Check if parameter set has given parameter


    .. cpp:function:: bool has_parameter_set(std::string key) const
    
        Check if parameter set has given nested parameter set


    .. cpp:function:: void get_parameter_keys(std::vector<std::string>& keys) const
    
        Return a vector of parameter keys


    .. cpp:function:: void get_parameter_set_keys(std::vector<std::string>& keys) const
    
        Return a vector of parameter set keys


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void parse_common(int argc, char* argv[])
    
        Parse filtered options (everything except PETSc options)


    .. cpp:function:: void parse_petsc(int argc, char* argv[])
    
        Parse filtered options (only PETSc options)


