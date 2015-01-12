
.. Documentation for the header file dolfin/plot/plot.h

.. _programmers_reference_cpp_plot_plot:

plot.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void interactive(bool really=false)

    Make the current plots interactive. If really is set, the interactive
    mode is entered even if 'Q' has been pressed.


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(const Variable&, std::string title="", std::string mode="auto")

    Plot variable of any supported type


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(std::shared_ptr<const Variable>, std::string title="", std::string mode="auto")

    Plot variable (shared_ptr version)


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(const Variable&, const Parameters& parameters)

    Plot variable (parameter version)


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(std::shared_ptr<const Variable>, std::shared_ptr<const Parameters> parameters)

    Plot variable (parameter, shared_ptr version)


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(const Expression& expression, const Mesh& mesh, std::string title="", std::string mode="auto")

    Plot expression


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(std::shared_ptr<const Expression> expression, std::shared_ptr<const Mesh> mesh, std::string title="", std::string mode="auto")

    Plot expression (shared_ptr version)


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(const Expression& expression, const Mesh& mesh, const Parameters& parameters)

    Plot expression (parameter version)


.. cpp:function:: std::shared_ptr<VTKPlotter> plot(std::shared_ptr<const Expression> expression, std::shared_ptr<const Mesh> mesh, std::shared_ptr<const Parameters> parameters)

    Plot expression (parameter, shared_ptr version)


