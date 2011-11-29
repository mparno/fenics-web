
Collection of documented demos
==============================

.. toctree::
   :glob:
   :numbered:
   :maxdepth: 1

   */*/cpp/documentation


To run the C++ demos, follow the below procedure:

* Download the source files i.e., ``main.cpp`` and ``CMakeLists.txt``,
  from the demo that you want to run. Some demos also contain UFL form
  files, e.g., ``Poisson.ufl``. Note that there may be multiple form
  files.

* Compile the form files to generate code with DOLFIN wrappers::

      $ ffc -l dolfin Poisson.ufl

  Generated .h files are usually distributed with the demos so you may
  choose to skip this step and use the provided header file directly,
  in this case ``Poisson.h``.

  If you wish to use optimized generated code, do::

      $ ffc -O -l dolfin Poisson.ufl

* Configure the demo build process::

      $ cmake .

* Compile the demo::

      $ make

* Run the created executable::

    $ ./demo



.. note::

    You must have a working installation of FEniCS in order to run the
    demos.

