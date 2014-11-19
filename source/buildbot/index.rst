

########
Buildbot
########

**************
About Buildbot
**************

The `Buildbot <http://www.buildbot.net>`_ is a system to automate the
compile/test cycle required by most software projects to validate code
changes. By automatically rebuilding and testing the tree each time
something has changed, build problems are pinpointed quickly, before
other developers are inconvenienced by the failure. The guilty developer
can be identified and harassed without human intervention. By running
the builds on a variety of platforms, developers who do not have the
facilities to test their changes everywhere before checkin will at least
know shortly afterwards whether they have broken the build or not.

***************
FEniCS buildbot
***************

The FEniCS buildbot is available from http://fenicsproject.org:8010. The
default `Waterfall <http://fenicsproject.org:8010/waterfall>`__ is quite
complex but the `Console View <http://fenicsproject.org:8010/console>`__
provides a good overview of the current status. In addition, the
following links might be useful when only some of the components are of
interest.

.. tabularcolumns:: |c|c|c|

.. list-table::
    :widths: 10, 10, 10
    :header-rows: 0
    :class: center

    * - Component
      - Waterfall
      - Console view
    * - :ref:`All components <about_components>`
      - `master <http://fenicsproject.org:8010/waterfall?category=dolfin.master&category=ferari.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=uflacs.master>`__
        / `next <http://fenicsproject.org:8010/waterfall?category=dolfin.next&category=ferari.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next&category=uflacs.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?category=dolfin.maint&category=ferari.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint&category=uflacs.maint>`__
        / `all <http://fenicsproject.org:8010/waterfall?category=dolfin.master&category=ferari.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=uflacs.master&category=dolfin.next&category=ferari.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next&category=uflacs.next&category=dolfin.maint&category=ferari.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint&category=uflacs.maint>`__
      - `master <http://fenicsproject.org:8010/console?category=dolfin.master&category=ferari.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=uflacs.master>`__
	/ `next <http://fenicsproject.org:8010/console?category=dolfin.next&category=ferari.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next&category=uflacs.next>`__
        / `maint <http://fenicsproject.org:8010/console?category=dolfin.maint&category=ferari.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint&category=uflacs.maint>`__
        / `all <http://fenicsproject.org:8010/console?category=dolfin.master&category=ferari.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=uflacs.master&category=dolfin.next&category=ferari.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next&category=uflacs.next&category=dolfin.maint&category=ferari.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint&category=uflacs.maint>`__
    * - :ref:`Core components <about_components_core>`
      - `master <http://fenicsproject.org:8010/waterfall?category=dolfin.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?category=dolfin.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next>`__
        / `maint <http://fenicsproject.org:8010/waterfall?category=dolfin.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint>`__
        / `all <http://fenicsproject.org:8010/waterfall?category=dolfin.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=dolfin.next&category=ffc.next&category=instant.next&category=ufc.next&category=ufl.next&category=dolfin.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint>`__
      - `master <http://fenicsproject.org:8010/console?category=dolfin.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master>`__
	/ `next <http://fenicsproject.org:8010/console?category=dolfin.next&category=ffc.next&category=fiat.next&category=instant.next&category=ufc.next&category=ufl.next>`__
        / `maint <http://fenicsproject.org:8010/console?category=dolfin.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint>`__
        / `all <http://fenicsproject.org:8010/console?category=dolfin.master&category=ffc.master&category=fiat.master&category=instant.master&category=ufc.master&category=ufl.master&category=dolfin.next&category=ffc.next&category=instant.next&category=ufc.next&category=ufl.next&category=dolfin.maint&category=ffc.maint&category=fiat.maint&category=instant.maint&category=ufc.maint&category=ufl.maint>`__

    * - :ref:`DOLFIN <about_components_dolfin>`
      - `master <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.master&category=dolfin.next&category=dolfin.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.master&category=dolfin.next&category=dolfin.maint>`__
    * - :ref:`FFC <about_components_ffc>`
      - `master <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.master&category=ffc.next&category=ffc.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=ffc&category=ffc.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=ffc&category=ffc.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=ffc&category=ffc.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=ffc&category=ffc.master&category=ffc.next&category=ffc.maint>`__
    * - :ref:`FIAT <about_components_fiat>`
      - `master <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.master&category=fiat.next&category=fiat.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=fiat&category=fiat.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=fiat&category=fiat.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=fiat&category=fiat.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=fiat&category=fiat.master&category=fiat.next&category=fiat.maint>`__
    * - :ref:`Instant <about_components_instant>`
      - `master <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.master&category=instant.next&category=instant.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=instant&category=instant.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=instant&category=instant.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=instant&category=instant.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=instant&category=instant.master&category=instant.next&category=instant.maint>`__
    * - :ref:`UFC <about_components_ufc>`
      - `master <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.master&category=ufc.next&category=ufc.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=ufc&category=ufc.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=ufc&category=ufc.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=ufc&category=ufc.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=ufc&category=ufc.master&category=ufc.next&category=ufc.maint>`__
    * - :ref:`UFL <about_components_ufl>`
      - `master <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.master&category=ufl.next&category=ufl.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=ufl&category=ufl.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=ufl&category=ufl.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=ufl&category=ufl.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=ufl&category=ufl.master&category=ufl.next&category=ufl.maint>`__
    * - uflacs
      - `master <http://fenicsproject.org:8010/waterfall?project=uflacs&category=uflacs.master>`__
	/ `next <http://fenicsproject.org:8010/waterfall?project=uflacs&category=uflacs.next>`__
	/ `maint <http://fenicsproject.org:8010/waterfall?project=uflacs&category=uflacs.maint>`__
	/ `all <http://fenicsproject.org:8010/waterfall?project=uflacs&category=uflacs.master&category=uflacs.next&category=uflacs.maint>`__
      - `master <http://fenicsproject.org:8010/console?project=uflacs&category=uflacs.master>`__
	/ `next <http://fenicsproject.org:8010/console?project=uflacs&category=uflacs.next>`__
	/ `maint <http://fenicsproject.org:8010/console?project=uflacs&category=uflacs.maint>`__
	/ `all <http://fenicsproject.org:8010/console?project=uflacs&category=uflacs.master&category=uflacs.next&category=uflacs.maint>`__

There are also other alternative renderings of the results:

* Feeds: `Atom <http://fenicsproject.org:8010/atom>`__ and `RSS
  <http://fenicsproject.org:8010/rss>`__
* `one line per build
  <http://fenicsproject.org:8010/one_line_per_build>`__ and `one box per
  builder <http://fenicsproject.org:8010/builders>`__
* Configuration: `all build slaves
  <http://fenicsproject.org:8010/buildslaves>`__ and `buildbot version
  information <http://fenicsproject.org:8010/about>`__


**********
Try server
**********

The Buildbot has a facility called "try server". This allows you to run
a build with your local changes before they are committed. To set up
your machine to run try builds, first install a recent version of
buildbot. Then go to your git-based working copy that contains changes
and run the following command::

    git diff | buildbot --patchlevel=1 \
                        --connect=pb \
                        --master=fenicsproject.org:8031 \
                        --username=<username> \
                        --password=<password> \
                        --who=<your name> \
                        --builder=<builder-name> \
                        --diff=-

To save some typing, add a file ``~/.buildbot/options`` with the following
contents::

    try_connect = 'pb'
    try_master = 'fenicsproject.org:8031'
    try_username = 'username'
    try_password = 'password'
    try_who = 'your name'

You can then start a build simply by running::

    git diff | buildbot try --patchlevel=1 --builder=<builder-name> --diff=-

To see a list of available options, see ``buildbot try --help``. For
instance, using ``--dryrun`` will gather info but not submit, while
using ``--get-builder-names`` will list the names of the available
builders that can be used with the ``--builder`` option. The builders
can also be set in ``~/.buildbot/options``, for instance::

    try_builders = ["dolfin-master-full-lucid-amd64", "dolfin-master-full-osx-10.7"]

For more information on running try builds, see the `try section
<http://buildbot.net/buildbot/docs/current/manual/cmdline.html#cmdline-try>`__
in the Buildbot documentation.

.. note::

    To be able to run try builds, you will need a username and
    password. This can be obtained by contacting `Johannes Ring
    <https://launchpad.net/~johannr>`__. For now, the access is limited
    to currently active developers.

.. warning::

    Doing try builds on the FEniCS buildbot is currently experimental
    and might not always works as expected.
