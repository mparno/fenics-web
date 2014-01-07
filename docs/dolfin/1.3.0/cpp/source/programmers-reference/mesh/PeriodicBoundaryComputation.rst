
.. Documentation for the header file dolfin/mesh/PeriodicBoundaryComputation.h

.. _programmers_reference_cpp_mesh_periodicboundarycomputation:

PeriodicBoundaryComputation.h
=============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PeriodicBoundaryComputation

    This class computes map from slave entity to master entity


    .. cpp:function:: static std::map<unsigned int, std::pair<unsigned int, unsigned int> > compute_periodic_pairs(const Mesh& mesh, const SubDomain& sub_domain, const std::size_t dim)
    
        For entities of dimension dim, compute map from a slave entity
        on this process (local index) to its master entity (owning
        process, local index on owner). If a master entity is shared
        by processes, only one of the owning processes is returned.


    .. cpp:function:: static MeshFunction<std::size_t> masters_slaves(boost::shared_ptr<const Mesh> mesh, const SubDomain& sub_domain, const std::size_t dim)
    
        This function returns a MeshFunction which marks mesh entities
        of dimension dim according to:
        
            2: slave entities
            1: master entities
            0: all other entities
        
        It is useful for visualising and debugging the Expression::map
        function that is used to apply periodic boundary conditions.


