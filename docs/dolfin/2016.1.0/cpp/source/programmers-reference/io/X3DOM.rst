
.. Documentation for the header file dolfin/io/X3DOM.h

.. _programmers_reference_cpp_io_x3dom:

X3DOM.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: X3DOMParameters

    Class data to store X3DOM view parameters.


    .. cpp:function:: enum class Representation
    
        X3DOM representation type


    .. cpp:function:: X3DOMParameters()
    
        Constructor (with default parameter settings)


    .. cpp:function:: void set_representation(Representation representation)
    
        Set representation of object (wireframe, surface or
        surface_with_edges)


    .. cpp:function:: Representation get_representation() const
    
        Get the current representation of the object (wireframe, surface or
        surface_with_edges)


    .. cpp:function:: std::array<double, 2> get_viewport_size() const
    
        Get the size of the viewport


    .. cpp:function:: void set_diffuse_color(std::array<double, 3> rgb)
    
        Set the RGB color of the object


    .. cpp:function:: std::array<double, 3> get_diffuse_color() const
    
        Get the RGB diffuse color of the object


    .. cpp:function:: void set_emissive_color(std::array<double, 3> rgb)
    
        Set the RGB emissive color


    .. cpp:function:: std::array<double, 3> get_emissive_color() const
    
        Get the RGB emissive color


    .. cpp:function:: void set_specular_color(std::array<double, 3> rgb)
    
        Set the RGB specular color


    .. cpp:function:: void set_background_color(std::array<double, 3> rgb)
    
        Set background RGB color


    .. cpp:function:: std::array<double, 3> get_background_color() const
    
        Get background RGB color


    .. cpp:function:: void set_ambient_intensity(double intensity)
    
        Set the ambient lighting intensity


    .. cpp:function:: double get_ambient_intensity() const
    
        Get the ambient lighting intensity


    .. cpp:function:: void set_shininess(double shininess)
    
        Set the surface shininess of the object


    .. cpp:function:: double get_shininess() const
    
        Set the surface shininess of the object


    .. cpp:function:: void set_transparency(double transparency)
    
        Set the transparency (0-1)


    .. cpp:function:: double get_transparency() const
    
        Get the transparency (0-1)


    .. cpp:function:: void set_color_map(const std::vector<double>& color_data)
    
        Set the color map by supplying a vector of 768 values
        (256*RGB) (using std::vector for Python compatibility via
        SWIG)


    .. cpp:function:: std::vector<double> get_color_map() const
    
        Get the color map as a vector of 768 values (256*RGB) (using
        std::vector for Python compatibility via SWIG)


    .. cpp:function:: boost::multi_array<float, 2> get_color_map_array() const
    
        Get the color map as a boost::multi_array (256x3)


    .. cpp:function:: void set_viewpoint_buttons(bool show)
    
        Toggle viewpoint buttons


    .. cpp:function:: bool get_viewpoint_buttons() const
    
        Get the viewpoint button state


    .. cpp:function:: void set_x3d_stats(bool show)
    
        Turn X3D 'statistics' window on/off


.. cpp:class:: X3DOM

    This class implements output of meshes to X3DOM XML or HTML5
    with X3DOM strings. The latter can be used for interactive visualisation
    
    Developer note: pugixml is used to created X3DOM and HTML5. By
    using pugixml, we produce valid XML, but care must be taken that
    the XML is also valid HTML. This includes not letting pugixml
    create self-closing elements, in cases. E.g., <foo
    bar="foobar"></foo> is fine, but the self-closing syntax <foo
    bar="foobar" /> while being valid XML is is not valid HTML5. See
    https://github.com/x3dom/x3dom/issues/600.


    .. cpp:function:: static std::string str(const Mesh& mesh, X3DOMParameters parameters=X3DOMParameters())
    
        Return X3D string for a Mesh


    .. cpp:function:: static std::string str(const Function& u, X3DOMParameters parameters=X3DOMParameters())
    
        Return X3D string for a Function


    .. cpp:function:: static std::string html(const Mesh& mesh, X3DOMParameters parameters=X3DOMParameters())
    
        Return HTML5 string with embedded X3D for a Mesh


    .. cpp:function:: static std::string html(const Function& u, X3DOMParameters parameters=X3DOMParameters())
    
        Return HTML5 string with embedded X3D for a Function


