.. General notes on how to run the C++ demos.

.. _demos_cpp_running_demos:

*****************
Running the demos
*****************

To run the ``C++`` demos, follow the below procedure:

* download the source files i.e., ``main.cpp`` and ``SConstruct``, from the
  demo that you want to run. Some demos also contain UFL form files e.g.,
  ``Poisson.ufl``. Note that there may be multiple form files.

* compile the form files files using DOLFIN as output language::

      $ ffc -l dolfin Poisson.ufl

* compile the ``main.cpp`` file against DOLFIN by running scons::

      $ scons

* run the created executable::

    $ ./demo

.. note::

    You must have a working installation of FEniCS in order to run the demo,
    see :ref:`installation` for details.

