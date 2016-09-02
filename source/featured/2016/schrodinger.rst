##########################################
Solving the Schrödinger Equation in FEniCS
##########################################

| Featured article 2016-09-02
| *Created by Petter Rosander, Anna Samuelsson and Nicklas Österbacka*

This article is a summary of a bachelor thesis conducted at Chalmers
University of Technology and the University of Gothenburg. The task at
hand is to solve the Schrödinger equation in FEniCS. In this
particular project the equation is applied to the physical system of a
hydrogen atom. Regarding time independent systems, the Schrödinger
equation is formulated as in [1].

.. math::

    \mathcal{\hat H} \psi(x) = E \psi(x),

whilst in the time dependent case, it is expressed as [2]

.. math::

   i \hbar \frac{\partial}{\partial t} \psi(x, t) = \mathcal{\hat H} \psi(x, t)

In the equations above, the eigenvalue $E$ is the systems energy,
:math:`\hbar` is Planck's reduced constant and the Hamiltonian
:math:`\mathcal{\hat H}` is defined as the sum of the systems kinetic,
as well as its potential energy. A normalization of :math:`|\psi|^2`
is a probability density function describing the location of a certain
particle [3]. In the case of the hydrogen atom, :math:`|\psi|^2`
defines the probability of finding the hydrogen electrone at a certain
position. :math:`\psi_{n,l,m}` denotes :math:`\psi` for a certain
quantum state where n is the principal quantum number, l is the
azimuthal quantum number and m is the magnetic quantum number
[3]. Certain linear cominations of :math:`\psi_{n,l,m}` are of certain
interest. These combinations are denoted :math:`\psi_{n,l}^w`.  The
potential energy of a hydrogen atom has a singularity in origo, which
is why a spherical scale geometry around origo is implemented in
FEniCS. Throughout this article, the distance from origo to the
spherical scale is referred to as the hole radius. The time dependence
has been approximated by Euler's backward method.

According to the results of this present project, the geometry of a
spherical scale is preferable to that of a sphere, as the smaller hole
radius tends to yield a slow and possibly false eigenvalue
convergence. Moreover, the numerically computed :math:`|\psi|^2` seems
to share the behaviour of its analytical equivalent, albeit incoherent
amplitude values. This could be solved by an increase in CPU
capacity. As Euler's backward method is A-stable [4], the numerically
computed :math:`\psi` approaches zero with time, even though its
analytical equivalent oscillates with a constant
amplitude[5]. Nevertheless, with time steps proportionate to the
length of the time interval, a concordant simulation can be obtained.

********************************
Simulations of the Hydrogen Atom
********************************

.. raw:: html

   <div class="container-fluid">
    <div class="row">
      <img src="../../_images/eig.png" class="img-responsive" alt="eig">
    </div>
  </div>

**Figure 1**: Numerically calculated eigenvalues of
:math:`\psi_{1,0,0}` in relation to mesh size. What differs the
calculations apart is the hole radiuses of the geometries. The red
line shows the analytically calculated eigenvalue.

.. raw:: html

   <div class="container-fluid">
    <div class="row">
      <img src="../../_images/eig2.png" class="img-responsive" alt="eig2">
    </div>
  </div>

**Figure 2**: To the left, an analytically calculated iso surfaces of
:math:`|\psi_{3,1}^w|^2` is found. The right picture shows its
numerically calculated equivalence.

.. raw:: html

   <div class="container-fluid">
    <div class="row">
      <img src="../../_images/eig3.png" class="img-responsive" alt="eig3">
    </div>
  </div>

**Figure 3**: Time developement of :math:`\psi_{1,0,0}`, numerically
calculated with different time step sizes.

For more information, see
http://publications.lib.chalmers.se/records/fulltext/238519/238519.pdf

******************************
References and acknowledgments
******************************

[1] Attila Szabo and Neil S. Ostlund. Modern Quantum Chemistry -
Introduction to Advanced Electronic Structure Theory. Dover
Publications, Inc., New York, 1st edition, 1996.

[2] Gunnar Ohlén. Kvantvärldens fenomen - teori och
begrepp. Studentlitteratur AB, Lund, 1st edition, 2005.

[3] Peter Atkins and Loretta Jones. Chemical Principles - The Quest
for Insight. W. H. Freeman and Company, New York, 5
edition, 2010.

[4] T. Lakoba. Stability analysis of finite-difference methods for
odes. Retrieved 2016- 04-14 from
http://www.cems.uvm.edu/tlakoba/math337/notes_4.pdf, p.45.

[5] Måns Henningson. Börja med kvantfysik. Unpublished book for the
course FUF040, Chalmers University of Technology, p. 79-80, 2015.
