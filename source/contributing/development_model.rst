.. _development_model:

#################
Development model
#################

***************
The big picture
***************

FEniCS development follows the model illustrated below.

.. image:: images/development_model.png
    :align: center

We note the following important points:

* Development happens in the **trunk**.
* At any time, there may exist one or more **stable branches**.
* **Stable releases** are made from stable branches.
* Blueprints may be targeted to a future **series** (for which there
  is yet no branch of the source repository).

***************
Stable releases
***************

The following process should be followed for the creation of a stable
release (a version X.Y.0). The primary requirement for making a
release is that the buildbot is green.

#. A **suggestion** is made on the FEniCS mailing list to create a new
   stable release.
#. Concensus is reached that the timing is good.
#. A **release manager** is appointed in consensus among the
   developers. The release manager is responsible for making sure that
   the steps below are carried out for each FEniCS component.
#. Create a Launchpad **series** X.Y.x and create a **milestone** for the
   release in the series.
#. Create a **branch** X.Y.x and associate it with the series.
#. Ensure that the buildbots are green.
#. For each project, run the fenics-release script and follow the
   instructions. The fenics-release script is available from

      bzr branch lp:fenics

   Note in particular the post-release instructions issued by the
   script.

A reasonable schedule for the release process is one week from the
time consensus is reached (2) to the time when the release is made
(7).

A release of X.Y.0 may be preceeded by X.Y-betaN and X.Y-rcN
releases. Once X.Y.0 has been released, it may be followed by bug fix
releases X.Y.1, X.Y.2, etc.

..
   *****************
   Snapshot releases
   *****************

   The creation of a stable release involves a fair amount of
   administration and it is also a lengthy process. Another type of
   release is a **snapshot release** made directly from trunk in between
   stable releases. Such an *ad hoc* release can be made at any time, as
   long as the buildbots are green.

   Developers may want to create snapshot releases for many reasons:

   * to point a collaborator to a fixed snapshot of trunk;
   * to get more testing of a new feature;
   * if the ChangeLog is growing long;
   * for the fun of making a release.

   The following simple procedure should be followed for snapshot
   releases:

   1. Announce the intention to make a snapshot release on the FEniCS mailing list.
   2. Wait a day.
   3. Make the release.

   Snapshot releases don't need to wait for new features to be completed
   as long as the buildbots are green. We can always make a new snapshot
   release when that feature has been implemented.

********************
Miscellaneous issues
********************

* A bug should only be marked as "fix committed" when it has been
  merged into trunk or a stable branch (not when it is pushed to a
  personal branch).
* X.Y-betaN releases are used to grind out bugs.
* X.Y-rcN releases are made when we believe there are "no" bugs.


