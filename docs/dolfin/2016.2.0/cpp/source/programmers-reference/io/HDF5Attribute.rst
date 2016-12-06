
.. Documentation for the header file dolfin/io/HDF5Attribute.h

.. _programmers_reference_cpp_io_hdf5attribute:

HDF5Attribute.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: HDF5Attribute

    HDF5Attribute gives access to the attributes of a dataset
    via set() and get() methods


    .. cpp:function:: HDF5Attribute(const hid_t hdf5_file_id, std::string dataset_name)
    
        Constructor


    .. cpp:function:: bool exists(const std::string attribute_name) const
    
        Check for the existence of an attribute on a dataset


    .. cpp:function:: void set(const std::string attribute_name, const double value)
    
        Set the value of a double attribute in the HDF5 file


    .. cpp:function:: void set(const std::string attribute_name, const std::size_t value)
    
        Set the value of a double attribute in the HDF5 file


    .. cpp:function:: void set(const std::string attribute_name, const std::vector<double>& value)
    
        Set the value of an array of float attribute in the HDF5 file


    .. cpp:function:: void set(const std::string attribute_name, const std::vector<std::size_t>& value)
    
        Set the value of an array of float attribute in the HDF5 file


    .. cpp:function:: void set(const std::string attribute_name, const std::string value)
    
        Set the value of a string attribute in the HDF5 file


    .. cpp:function:: void get(const std::string attribute_name, double& value) const
    
        Set the value of a double attribute in the HDF5 file


    .. cpp:function:: void get(const std::string attribute_name, std::vector<double>& value) const
    
        Get the value of a vector double attribute in the HDF5 file


    .. cpp:function:: void get(const std::string attribute_name, std::size_t& value) const
    
        Set the value of a double attribute in the HDF5 file


    .. cpp:function:: void get(const std::string attribute_name, std::vector<std::size_t>& value) const
    
        Get the value of a vector double attribute in the HDF5 file


    .. cpp:function:: void get(const std::string attribute_name, std::string& value) const
    
        Get the value of an attribute in the HDF5 file as a string


    .. cpp:function:: const std::string str(const std::string attribute_name) const
    
        Get the value of the attribute in the HDF5 file
        as a string representation


    .. cpp:function:: const std::string type_str(const std::string attribute_name) const
    
        Get the type of the attribute "string", "float", "int"
        "vectorfloat", "vectorint" or "unsupported"


    .. cpp:function:: const std::string str() const
    
        Get the names of all the attributes on this dataset


    .. cpp:function:: const std::vector<std::string> list_attributes() const
    
        Get the names of all the attributes on this dataset as a
        std::vector<std::string>


