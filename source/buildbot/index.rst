

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
      - `1.0.x <http://fenicsproject.org:8010/waterfall?category=dolfin.1.0.x&category=ferari.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=syfi.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x&category=viper.1.0.x>`__
        / `trunk <http://fenicsproject.org:8010/waterfall?category=dolfin.trunk&category=ferari.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=syfi.trunk&category=ufc.trunk&category=ufl.trunk>`__
        / `both <http://fenicsproject.org:8010/waterfall?category=dolfin.1.0.x&category=ferari.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=syfi.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x&category=viper.1.0.x&category=dolfin.trunk&category=ferari.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=syfi.trunk&category=ufc.trunk&category=ufl.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?category=dolfin.1.0.x&category=ferari.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=syfi.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x&category=viper.1.0.x>`__
        / `trunk <http://fenicsproject.org:8010/console?category=dolfin.trunk&category=ferari.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=syfi.trunk&category=ufc.trunk&category=ufl.trunk>`__
        / `both <http://fenicsproject.org:8010/console?category=dolfin.1.0.x&category=ferari.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=syfi.trunk&category=ufc.2.0.x&category=ufl.1.0.x&category=viper.trunk&category=dolfin.trunk&category=ferari.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=syfi.trunk&category=ufc.trunk&category=ufl.trunk>`__
    * - :ref:`Core components <about_components_core>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?category=dolfin.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x>`__
        / `trunk <http://fenicsproject.org:8010/waterfall?category=dolfin.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=ufc.trunk&category=ufl.trunk>`__
        / `both <http://fenicsproject.org:8010/waterfall?category=dolfin.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x&category=dolfin.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=ufc.trunk&category=ufl.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?category=dolfin.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x>`__
        / `trunk <http://fenicsproject.org:8010/console?category=dolfin.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=ufc.trunk&category=ufl.trunk>`__
        / `both <http://fenicsproject.org:8010/console?category=dolfin.1.0.x&category=ffc.1.0.x&category=fiat.1.0.x&category=instant.1.0.x&category=ufc.2.0.x&category=ufl.1.0.x&category=dolfin.trunk&category=ffc.trunk&category=fiat.trunk&category=instant.trunk&category=ufc.trunk&category=ufl.trunk>`__

    * - :ref:`DOLFIN <about_components_dolfin>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=dolfin&category=dolfin.1.0.x&category=dolfin.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=dolfin&category=dolfin.1.0.x&category=dolfin.trunk>`__
    * - :ref:`FFC <about_components_ffc>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=ffc&category=ffc.1.0.x&category=ffc.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=ffc&category=ffc.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=ffc&category=ffc.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=ffc&category=ffc.1.0.x&category=ffc.trunk>`__
    * - :ref:`FIAT <about_components_fiat>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=fiat&category=fiat.1.0.x&category=fiat.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=fiat&category=fiat.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=fiat&category=fiat.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=fiat&category=fiat.1.0.x&category=fiat.trunk>`__
    * - :ref:`Instant <about_components_instant>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=instant&category=instant.1.0.x&category=instant.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=instant&category=instant.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=instant&category=instant.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=instant&category=instant.1.0.x&category=instant.trunk>`__
    * - :ref:`SyFi <about_components_syfi>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=syfi&category=syfi.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=syfi&category=syfi.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=syfi&category=syfi.1.0.x&category=syfi.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=syfi&category=syfi.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=syfi&category=syfi.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=syfi&category=syfi.1.0.x&category=syfi.trunk>`__
    * - :ref:`UFC <about_components_ufc>`
      - `2.0.x <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.2.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=ufc&category=ufc.2.0.x&category=ufc.trunk>`__
      - `2.0.x <http://fenicsproject.org:8010/console?project=ufc&category=ufc.2.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=ufc&category=ufc.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=ufc&category=ufc.2.0.x&category=ufc.trunk>`__
    * - :ref:`UFL <about_components_ufl>`
      - `1.0.x <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.trunk>`__
	/ `both <http://fenicsproject.org:8010/waterfall?project=ufl&category=ufl.1.0.x&category=ufl.trunk>`__
      - `1.0.x <http://fenicsproject.org:8010/console?project=ufl&category=ufl.1.0.x>`__
	/ `trunk <http://fenicsproject.org:8010/console?project=ufl&category=ufl.trunk>`__
	/ `both <http://fenicsproject.org:8010/console?project=ufl&category=ufl.1.0.x&category=ufl.trunk>`__

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
buildbot. Then go to your bzr-based working copy that contains changes
and run the following command::

    bzr diff | buildbot --connect=pb \
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

    bzr diff | buildbot try --builder=<builder-name> --diff=-

To see a list of available options, see ``buildbot try --help``. For
instance, using ``--dryrun`` will gather info but not submit, while
using ``--get-builder-names`` will list the names of the available
builders that can be used with the ``--builder`` option. The builders
can also be set in ``~/.buildbot/options``, for instance::

    try_builders = ["dolfin-trunk-full-lucid-amd64", "dolfin-trunk-full-osx-10.7"]

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
