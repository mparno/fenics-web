
##########################################################################
Subject-specific simulations of cerebrospinal fluid flow and drug delivery
##########################################################################

| Featured article 2015-06-12
| *Created by P.T Haga, M. Mortensen and K.-A. Mardal*

The cerebrospinal fluid (CSF) surrounds the central nervous system (CNS), and drugs infused into the CSF can thus quickly be absorbed. However, because of large subject-specific variations and a complex oscillating flow, it has been proven difficult to predict and control the addition of drugs [1]_. A study has been done to simulate the flow and drug transport of the CSF in a subject-specific model of the cranial region of the spine using FEniCS. This article briefly describes the results of the study.

***********************************
Subject-specific computational mesh
***********************************
A three dimensional geometry of the cervical subarachnoid space is shown in Figure 1. The geometry was provided by Dr. Bryn Martin at the Conquer Chiari Research Center at the University of Akron, and more information can be found in [2]_. The anatomical model was constructed based on manual segmentation of T2-weighted magnetic resonance (MR) image sequences of a healthy volunteer using freely available software ITK-Snap (Version 2.2, University of Pennsylvania). Idealized NRDL were separately constructed and added to the model using Autodesk Maya (Autodesk Inc., Cleveland, OH). Nonuniform unstructured computational meshes were generated using ANSYS ICEM CFD (ANSYS Inc., Canonsburg, PA). An example of the surface elements of the computational mesh is shown in Figure 1 (c). The complete geometry [Figure 1 (a)] is 18 cm from top to bottom and the two end planes are both placed in the xz-plane.

.. image:: images/geometry.png
  :scale: 50 %

**Figure 1:** SSS geometry and computational surface mesh. (a) Complete surface model. (b) Transparent model showing nerve roots and denticulate ligaments. (c) Surface elements of the final computational mesh.

****************************
Lagrangian Particle Tracking
****************************
Due to very low diffusivity of the drugs, numerical issues may arise when using the Finite Element Method on the advection-diffusion equation. Lagrangian Particle Tracking is a method where such issues are not present. This method uses discrete particles and utilizes the velocity field to compute the position of the particles. In other words, for every particle, we solve the equation

.. math::

  	\frac{\partial x_p}{\partial t} = u(x_p,t),

where :math:`x_p` is the position of one particle, and :math:`u` is the velocity field. The particle density in one cell is computed simply by taking

.. math::

	\rho_{\kappa} = \frac{\text{number of particles in cell}}{\text{volume of cell}}.

***************
CFD simulations
***************
The simulations were done using the Open Source high-performance Navier-Stokes solver Oasis [3]_ coupled with the Lagrangian Particle Tracking for FEniCS [4]_ made by M. Kuchta and M. Mortensen. The solver was run on Abel supercomputer using 96 cores. The resulting velocity field revealed the formation of vortical structures in conjuction with the nerve roots and denticulate ligaments, as shown in Figure 2.

.. image:: images/streamlines_csf.png
	:align: center

**Figure 2:** Streamlines of the cerebrospinal fluid shows vorticity around the nerve roots and denticulate ligaments.

The particles were injected with a velocity coinciding with a 1 ml injection over 1 minute with a 22 gauge needle. 54 particles were injected every time-step giving a total of approximately 130.000 particles after 5 cardiac cycles. Figure 3 shows an animation of the drug concentration during the 5 first cardiac cycles. A sensitivity analysis revealed a relatively large difference on the drug distribution for different posterior injection points and inclination, and a relatively small difference on different lateral injection points, and different injection velocities. Longer simulations are needed to establish the long-term sensitivity of these parameters.  


.. image:: images/scalar_anim.gif
	:align: center

**Figure 3:** An animation of the drug concentration for the first 5 cardiac cycles. The color represents the number of particles per :math:`mm^2`.

References
*************************************************************************


.. [1] G. Hocking and J. A. W. Wildsmith. Intrathecal drug spread. British Journal of Anaesthesia, 93(4):568–578, 2004. doi: 10.1093/bja/aeh204. URL http://bja.oxfordjournals.org/content/93/4/568.short.

.. [2] Soroush Heidari Pahlavian, Theresia Yiallourou, R. Shane Tubbs, Alexander C. Bunck, Francis Loth, Mark Goodin, Mehrdad Raisee, and Bryn A. Martin. The impact of spinal cord nerve roots and denticulate ligaments on cerebrospinal fluid dynamics in the cervical spine. PLoS ONE, 9(4):e91888, 04 2014. doi: 10.1371/journal.pone.0091888. URL http://dx.doi.org/10.1371%2Fjournal.pone.0091888.

.. [3] Mikael Mortensen and Kristian Valen-Sendstad. Oasis: A high- level/high-performance open source navier–stokes solver. Computer Physics Communications, 188(0):177 – 188, 2015. ISSN 0010-4655. doi: http://dx.doi.org/10.1016/j.cpc.2014.10.026. URL http://www.sciencedirect.com/science/article/pii/S0010465514003786.

.. [4] https://github.com/MiroK/lagrangian-particles

