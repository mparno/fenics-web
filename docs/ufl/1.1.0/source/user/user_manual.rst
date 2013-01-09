.. UFL user manual

.. _ufl_user_manual:

###############
UFL user manual
###############

.. note::

   This manual was was copied from LaTeX. It needs
   updating and reformatting. (MER, 09.06.2011)

.. toctree::
   :maxdepth: 1

   introduction
   form_language
   examples
   internal_representation
   algorithms
   command_line_utils
   installation

*****************
About this manual
*****************

Intended audience
=================

This manual is written both for the beginning and the advanced user.
There is also some useful information for developers. More advanced
topics are treated at the end of the manual or in the appendix.

Typographic conventions
=======================

Code is written in monospace ``like this``. Commands that should be
entered in a Unix shell::

  # ./configure
  # make

Commands are written in the dialect of the ``bash`` shell. For other
shells, such as ``tcsh``, appropriate translations may be needed.



Enumeration and list indices
============================

Throughout this manual, elements :math:`x_i` of sets :math:`\{x_i\}`
of size :math:`n` are enumerated from :math:`i = 0` to :math:`i =
n-1`. Derivatives in :math:`\mathbb{R}^n` are enumerated similarly
:math:`\frac{\partial}{\partial x_0}, \frac{\partial}{\partial x_1}, \ldots, \frac{\partial}{\partial x_{n-1}}`.

Contact
=======

Comments, corrections and contributions to this manual are most welcome
and should be sent to ufl@lists.launchpad.net

