__authors__ = "Anders Logg <logg@simula.no>, Kristian Oelgaard <k.b.oelgaard@gmail.com> and Marie E. Rognes <meg@simula.no>"
__copyright__ = "Copyright (C) 2010- " + __authors__
__license__  = "GNU LGPL version 3 or any later version"

library_index_cpp_template = """\
.. _programmers_reference_index:

############################################
C++ Programmer's reference for DOLFIN-%s
############################################

.. toctree::
   :glob:
   :maxdepth: 2

   */index
"""

module_index_cpp_template = """\
.. _programmers_reference_%s_index:

**********************************
DOLFIN %s module
**********************************

.. toctree::
   :glob:
   :maxdepth: 1

   *
"""

header_cpp_template = """\

.. Documentation for the header file dolfin/%(module)s/%(header)s

.. _programmers_reference_cpp_%(module)s_%(prefix)s:

%(title)s
%(formatting)s

.. note::
%(note)s

%(text)s
"""

note_template = """
The documentation on this page was automatically extracted from the
DOLFIN C++ code and may need to be edited or expanded.
"""

automodule_template = """
.. automodule:: %s
   :members:
   :undoc-members:
   :show-inheritance:
"""

def generate_library_index(version):
    return library_index_cpp_template % version

def generate_module_index(module):
    return module_index_cpp_template % (module, module)


