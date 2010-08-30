.. Documentation for the header file dolfin/plot/plot.h

.. _programmers_reference_cpp_plot_plot:

plot.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: void plot(const Expression& v, const Mesh& mesh,
                                std::string title="Expression", std::string mode="auto")
    
        Plot function

    .. cpp:function:: void plot(const Function& v,
                       std::string title="Function", std::string mode="auto")
    
        Simple built-in plot commands for plotting functions and meshes.
        For plotting to work, PyDOLFIN and Viper must be installed.
        Plot function

    .. cpp:function:: void plot(const Mesh& mesh,
                                std::string title="Mesh")
    
        Plot mesh

    .. cpp:function:: void plot(const MeshFunction<bool>& f,
                                std::string title="MeshFunction<bool>")
    
        Plot mesh function

    .. cpp:function:: void plot(const MeshFunction<double>& f,
                                std::string title="MeshFunction<double>")
    
        Plot mesh function

    .. cpp:function:: void plot(const MeshFunction<uint>& f,
                                std::string title="DOLFIN MeshFunction<uint>")
    
        Plot mesh function

