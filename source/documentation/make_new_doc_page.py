page = """
.. title:: Documentation for %(version)s

.. _overview_of_documentation_for_FEniCS_%(version)s:

####################################################
Documentation for FEniCS %(version)s
####################################################

Our documentation includes a book(!), a collection of documented demo
programs, and complete references for the FEniCS application
programming interface (API). Note that the FEniCS API is documented
separately for each FEniCS component. The most important interfaces
are those of the C++/Python problem solving environment :ref:`DOLFIN
<about_components_dolfin>` and the form language :ref:`UFL
<about_components_ufl>`.

(*This page accesses the FEniCS %(version)s documentation. Not the
version you are looking for? See* :ref:`other_versions_%(version)s`.)

.. raw:: html

  <div id="main">
  <div id="container" class="feature">
    <div id="content">
      <div id="sub-feature">
	<div id="front-block-1" class="front-block block">
          <h2> The FEniCS Tutorial </h2>

A good starting point for new users is the :doc:`FEniCS Tutorial
<tutorial/index>`. The tutorial will help you get quickly up and
running with solving differential equations in FEniCS. The tutorial
focuses exclusively on the FEniCS Python interface, since this is the
simplest approach to exploring FEniCS for beginners.

.. raw:: html

    <h2>The FEniCS Book</h2>
    <a href="/book/"><img alt=''src='/_static/images/fenics_book_cover.png' class='avatar avatar-84 photo' width='84'/></a>

:ref:`The FEniCS Book <book>`, *Automated Solution of Differential
Equations by the Finite Element Method*, is a comprehensive (700
pages) book documenting the mathematical methodology behind the FEniCS
Project and the software developed as part of the FEniCS Project. The
FEniCS Tutorial is included as the opening chapter of the FEniCS Book.

.. raw:: html

  <h2>The FEniCS Manual</h2>

`The FEniCS Manual
<http://launchpad.net/fenics-book/trunk/final/+download/fenics-manual-2011-10-31.pdf>`__
is a 200-page excerpt from the FEniCS Book, including the FEniCS
Tutorial, an introduction to the finite element method and
documentation of DOLFIN and UFL.

.. raw:: html

  </div><!-- #front-block-1 .front-block .block-->

  <div id="front-block-2" class="front-block block">
    <h2>Demos</h2>

A simple way to build your first FEniCS application is to copy and
modify one of the existing demos:

* `Documented DOLFIN demos (Python) <../documentation/dolfin/%(version_dot)s/python/demo/index.html>`__
* `Documented DOLFIN demos (C++) <../documentation/dolfin/%(version_dot)s/cpp/demo/index.html>`__

More demos can be found either in the ``demo`` directory of the DOLFIN
source tree or under one of the directories ``/usr/share/dolfin/demo``
(GNU/Linux),
``/Applications/FEniCS.app/Contents/Resources/share/dolfin/demo`` (Mac
OS X) or ``C:\FEniCS\share\dolfin\demo`` (Windows) if FEniCS was
installed from a binary package.

.. raw:: html

  <h2>Complete Programmer's References</h2>

* `All classes and functions in DOLFIN (Python) <../documentation/dolfin/%(version_dot)s/python/genindex.html>`__
* `All classes and functions in DOLFIN (C++) <../documentation/dolfin/%(version_dot)s/cpp/genindex.html>`__
* `All classes and functions in UFL <../documentation/ufl/%(ufl_version)s/genindex.html>`__

.. raw:: html

  <h2>Quick Programmer's Reference</h2>

We are working on adding a quick reference for common classes and
functions. It will be available here *soon*.


.. raw:: html

   </div><!-- #front-block-2 .front-block .block-->
   </div><!-- #sub-feature -->
     </div><!-- #content -->
       </div><!-- #container .feature -->
         </div><!-- #main -->



*************
Release notes
*************

If you are updating your application code to a new FEniCS release,
make sure to check the :ref:`release notes <releases>` where you will
find detailed information about new features and interface changes.

.. _other_versions_%(version)s:

***********************************
All versions of the documentation
***********************************

%(all_versions)s

.. toctree::
   :hidden:

   tutorial/index

"""

def write_doc_page(options):

    filename = "doc_%s.rst" % options["version"]
    file = open(filename, 'w')
    print "Generating %s" % filename
    contents = page % options
    file.write(contents)
    file.close()

def generate_doc_pages(doc_versions):

    all_versions = "\n".join(["* :doc:`doc_%s`" % v[0] for v in doc_versions])

    for (dolfin, ufl) in doc_versions:
        versions = {"version": dolfin,
                    "version_dot": dolfin.replace("-", "."),
                    "ufl_version": ufl,
                    "all_versions": all_versions}
        write_doc_page(versions)

if __name__ == "__main__":

    # Generate main documentation page(s) by list of (dolfin-version,
    # ufl-version):
    versions = [("dev", "dev"),
                ("1.0-rc2", "1.0-rc1"),
                ("1.0-rc1", "1.0-rc1"),
                ("1.0-beta2", "1.0-beta3"),
                ("1.0-beta", "1.0-beta2")]

    # Or just the newest version(s):
    # versions = [("1.0.0", "1.0.0")]
    generate_doc_pages(versions)
