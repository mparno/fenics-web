
.. Documentation for the header file dolfin/parameter/Parameter.h

.. _programmers_reference_cpp_parameter_parameter:

Parameter.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Parameter

    Base class for parameters.


    .. cpp:function:: explicit Parameter(std::string key)
    
        Create parameter for given key


    .. cpp:function:: std::string key() const
    
        Return parameter key


    .. cpp:function:: std::string description() const
    
        Return parameter description


    .. cpp:function:: bool is_set() const
    
        Return true if parameter is set, return false otherwise


    .. cpp:function:: void reset()
    
        Reset the parameter to empty, so that is_set() returns false.


    .. cpp:function:: std::size_t access_count() const
    
        Return access count (number of times parameter has been accessed)


    .. cpp:function:: std::size_t change_count() const
    
        Return change count (number of times parameter has been changed)


    .. cpp:function:: void set_range(int min_value, int max_value)
    
        Set range for int-valued parameter


    .. cpp:function:: void set_range(double min_value, double max_value)
    
        Set range for double-valued parameter


    .. cpp:function:: void set_range(std::set<std::string> range)
    
        Set range for string-valued parameter


    .. cpp:function:: void get_range(int& min_value, int& max_value) const
    
        Get range for int-valued parameter


    .. cpp:function:: void get_range(double& min_value, double& max_value) const
    
        Get range for double-valued parameter


    .. cpp:function:: void get_range(std::set<std::string>& range) const
    
        Get range for string-valued parameter


    .. cpp:function:: const Parameter& operator= (int value)
    
        Assignment from int


    .. cpp:function:: const Parameter& operator= (double value)
    
        Assignment from double


    .. cpp:function:: const Parameter& operator= (std::string value)
    
        Assignment from string


    .. cpp:function:: const Parameter& operator= (const char* value)
    
        Assignment from string


    .. cpp:function:: const Parameter& operator= (bool value)
    
        Assignment from bool


    .. cpp:function:: operator int() const
    
        Cast parameter to int


    .. cpp:function:: operator std::size_t() const
    
        Cast parameter to std::size_t


    .. cpp:function:: operator double() const
    
        Cast parameter to double


    .. cpp:function:: operator std::string() const
    
        Cast parameter to string


    .. cpp:function:: operator bool() const
    
        Cast parameter to bool


    .. cpp:function:: std::string type_str() const = 0
    
        Return value type string


    .. cpp:function:: std::string value_str() const = 0
    
        Return value string


    .. cpp:function:: std::string range_str() const = 0
    
        Return range string


    .. cpp:function:: std::string str() const = 0
    
        Return short string description


.. cpp:class:: IntParameter

    *Parent class(es)*
    
        * :cpp:class:`Parameter`
        
    Parameter with value type int


    .. cpp:function:: explicit IntParameter(std::string key)
    
        Create unset int-valued


    .. cpp:function:: IntParameter(std::string key, int value)
    
        Create int-valued parameter


    .. cpp:function:: void set_range(int min_value, int max_value)
    
        Set range


    .. cpp:function:: void get_range(int &min_value, int &max_value) const
    
        Get range


    .. cpp:function:: const IntParameter& operator= (int value)
    
        Assignment


    .. cpp:function:: operator int() const
    
        Cast parameter to int


    .. cpp:function:: operator std::size_t() const
    
        Cast parameter to std::size_t


    .. cpp:function:: std::string type_str() const
    
        Return value type string


    .. cpp:function:: std::string value_str() const
    
        Return value string


    .. cpp:function:: std::string range_str() const
    
        Return range string


    .. cpp:function:: std::string str() const
    
        Return short string description


.. cpp:class:: DoubleParameter

    *Parent class(es)*
    
        * :cpp:class:`Parameter`
        
    Parameter with value type double


    .. cpp:function:: explicit DoubleParameter(std::string key)
    
        Create unset double-valued parameter


    .. cpp:function:: DoubleParameter(std::string key, double value)
    
        Create double-valued parameter


    .. cpp:function:: void set_range(double min_value, double max_value)
    
        Set range


    .. cpp:function:: void get_range(double &min_value, double &max_value) const
    
        Get range


    .. cpp:function:: const DoubleParameter& operator= (double value)
    
        Assignment


    .. cpp:function:: operator double() const
    
        Cast parameter to double


    .. cpp:function:: std::string type_str() const
    
        Return value type string


    .. cpp:function:: std::string value_str() const
    
        Return value string


    .. cpp:function:: std::string range_str() const
    
        Return range string


    .. cpp:function:: std::string str() const
    
        Return short string description


.. cpp:class:: StringParameter

    *Parent class(es)*
    
        * :cpp:class:`Parameter`
        
    Parameter with value type string


    .. cpp:function:: explicit StringParameter(std::string key)
    
        Create unset string-valued parameter


    .. cpp:function:: StringParameter(std::string key, std::string value)
    
        Create string-valued parameter


    .. cpp:function:: void set_range(std::set<std::string> range)
    
        Set range


    .. cpp:function:: void get_range(std::set<std::string>& range) const
    
        Get range


    .. cpp:function:: const StringParameter& operator= (std::string value)
    
        Assignment


    .. cpp:function:: const StringParameter& operator= (const char* value)
    
        Assignment


    .. cpp:function:: operator std::string() const
    
        Cast parameter to string


    .. cpp:function:: std::string type_str() const
    
        Return value type string


    .. cpp:function:: std::string value_str() const
    
        Return value string


    .. cpp:function:: std::string range_str() const
    
        Return range string


    .. cpp:function:: std::string str() const
    
        Return short string description


.. cpp:class:: BoolParameter

    *Parent class(es)*
    
        * :cpp:class:`Parameter`
        
    Parameter with value type bool


    .. cpp:function:: explicit BoolParameter(std::string key)
    
        Create unset bool-valued parameter


    .. cpp:function:: BoolParameter(std::string key, bool value)
    
        Create bool-valued parameter


    .. cpp:function:: const BoolParameter& operator= (bool value)
    
        Assignment


    .. cpp:function:: operator bool() const
    
        Cast parameter to bool


    .. cpp:function:: std::string type_str() const
    
        Return value type string


    .. cpp:function:: std::string value_str() const
    
        Return value string


    .. cpp:function:: std::string range_str() const
    
        Return range string


    .. cpp:function:: std::string str() const
    
        Return short string description


