
.. Documentation for the header file dolfin/plot/VTKPlotter.h

.. _programmers_reference_cpp_plot_vtkplotter:

VTKPlotter.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: VTKPlotter

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class enables visualization of various DOLFIN entities.
    It supports visualization of meshes, functions, expressions, boundary
    conditions and mesh functions. It can plot data wrapped in classes
    conforming to the GenericVTKPlottable interface.
    The plotter has several parameters that the user can set and adjust to
    affect the appearance and behavior of the plot.
    
    A plotter can be created and used in the following way:
    
      Mesh mesh = ...;
      VTKPlotter plotter(mesh);
      plotter.plot();
    
    Parameters can be adjusted at any time and will take effect on the next
    call to the plot() method. The following parameters exist:
    
    ============== ============ ================ ====================================
     Name           Value type   Default value              Description
    ============== ============ ================ ====================================
     mode            String        "auto"         For vector valued functions,
                                                  this parameter may be set to
                                                  "glyphs" or "displacement".
                                                  Scalars may be set to "warp" in
                                                  2D only. A value of "color" is
                                                  valid in all cases; for vectors,
                                                  the norms are used. See below for
                                                  a summary of default modes,
                                                  used when set to "auto".
     interactive     Boolean     False            Enable/disable interactive mode
                                                  for the rendering window.
                                                  For repeated plots of the same
                                                  object (animated plots), this
                                                  parameter should be set to false.
     wireframe       Boolean     True for         Enable/disable wireframe
                                 meshes, else     rendering of the object.
                                 false
     title           String      Inherited        The title of the rendering
                                 from the         window
                                 name/label of
                                 the object
     scale           Double      1.0              Adjusts the scaling of the
                                                  warping and glyphs
     scalarbar       Boolean     False for        Hide/show the colormapping bar
                                 meshes, else
                                 true
     axes            Boolean     False            Show X-Y-Z axes.
    
     rescale         Boolean     True             Enable/disable recomputation
                                                  of the scalar to color mapping
                                                  on every iteration when performing
                                                  repeated/animated plots of the same
                                                  data. If both range_min and
                                                  range_max are set, this parameter
                                                  is ignored.
     range_min       Double                       Set lower range of data values.
                                                  Disables automatic (re-)computation
                                                  of the lower range.
     range_max       Double                       Set upper range of data values.
                                                  Disables automatic (re-)computation
                                                  of the upper range.
     elevate         Double      -65.0 for 2D     Set camera elevation.
                                 warped scalars,
                                 0.0 otherwise
     prefix          String      "dolfin_plot_"   Filename prefix used when
                                                  saving plots to file in
                                                  interactive mode. An integer
                                                  counter is appended after the
                                                  prefix.
     helptext        Boolean     True             Enable/disable the hover-over
                                                  help-text in interactive
                                                  mode
     window_width    Integer     600              The width of the plotting window
                                                  in pixels
     window_height   Integer     400              The height of the plotting window
                                                  in pixels
     tile_windows    Boolean     True             Automatically tile plot windows.
    
     key             String                       Key (id) of the plot window, used to
                                                  decide if a new plotter should be
                                                  created or a current one updated
                                                  when called through the static
                                                  plot() interface (in plot.h).
                                                  If not set, the object's unique
                                                  id (Variable::id) is used.
     input_keys      String      ""               Synthesize key presses, as if these
                                                  keys are pressed by the user in
                                                  the plot window.
                                                  For example: "ww++m" shows the data
                                                  as large points on a wireframe
                                                  mesh.
     hide_above      Double                       If either of these are set, scalar
     hide_below      Double                       values above or below will not be
                                                  shown in the plot.
    ============== ============ ================ ====================================
    
    The default visualization mode for the different plot types are as follows:
    
    =========================  ============================ =====================
     Plot type                  Default visualization mode   Alternatives
    =========================  ============================ =====================
     Meshes                     Wireframe rendering          None
     2D scalar functions        Scalar warping               Color mapping
     3D scalar functions        Color mapping                None
     2D/3D vector functions     Glyphs (vector arrows)       Displacements,
                                                             Color mapping (norm)
    =========================  ============================ =====================
    
    Expressions and boundary conditions are also visualized according to the
    above table.


    .. cpp:function:: VTKPlotter(boost::shared_ptr<const Variable>, QVTKWidget *widget = NULL)
    
        Create plotter for a variable. If a widget is supplied, this widget
        will be used for drawing, instead of a new top-level widget. Ownership
        is transferred.


    .. cpp:function:: VTKPlotter(boost::shared_ptr<const Expression> expression, boost::shared_ptr<const Mesh> mesh, QVTKWidget *wiget = NULL)
    
        Create plotter for an Expression with associated Mesh. If a widget is
        supplied, this widget will be used for drawing, instead of a new
        top-level widget. Ownership is transferred.


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void plot(boost::shared_ptr<const Variable> variable=boost::shared_ptr<const Variable>())
    
        Plot the object


    .. cpp:function:: void interactive(bool enter_eventloop = true)
    
        Make the current plot interactive


    .. cpp:function:: void write_png(std::string filename="")
    
        Save plot to PNG file (file suffix appended automatically, filename
        optionally built from prefix)


    .. cpp:function:: void write_pdf(std::string filename="")
    
        Save plot to PDF file (file suffix appended automatically, filename
        optionally built from prefix)


    .. cpp:function:: const std::string& key() const
    
        Return key (i.e., plotter id) of the object to plot


    .. cpp:function:: void set_key(std::string key)
    
        Set the key (plotter id)


    .. cpp:function:: static std::string to_key(const Variable &var)
    
        Return default key (plotter id) of a Variable (object to plot).


    .. cpp:function:: void azimuth(double angle)
    
        Camera control


    .. cpp:function:: static void all_interactive(bool really=false)
    
        Make all plot windows interactive. If really is set, the interactive
        mode is entered even if 'Q' has been pressed.


