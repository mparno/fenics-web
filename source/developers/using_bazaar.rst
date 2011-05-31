
.. _developers_using_bazaar:

************
Using Bazaar
************

Here is a quick reference for `using Bazaar
<http://doc.bazaar-vcs.org/bzr.2.0/en/quick-reference/index.html>`_.
In addition, a few useful commands for Bazaar follow below.

To set your identity with Bazaar, type

.. code-block:: sh

    bzr whoami "My Name <myname@foo.com>"

To create a new branch:

.. code-block:: sh

    bzr branch <address-to-original-branch> [address-to-new-branch]

To commit a change:

.. code-block:: sh

    bzr commit

To push changes:

.. code-block:: sh

    bzr push <address-to-branch>

To pull changes:

.. code-block:: sh

    bzr pull <address-to-branch>

The current development version of each FEniCS component can be
obtained directly using a special shortcut for code hosted on
Launchpad:

.. code-block:: sh

    bzr branch lp:<project-name>

For instance, one may create a branch of the main DOLFIN repository by
typing

.. code-block:: sh

    bzr branch lp:dolfin
