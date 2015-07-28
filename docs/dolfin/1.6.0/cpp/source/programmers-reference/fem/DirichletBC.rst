
.. Documentation for the header file dolfin/fem/DirichletBC.h

.. _programmers_reference_cpp_fem_dirichletbc:

DirichletBC.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DirichletBC

    *Parent class(es)*
    
        * :cpp:class:`Hierarchical<DirichletBC>`
        
        * :cpp:class:`Variable`
        
    This class specifies the interface for setting (strong)
    Dirichlet boundary conditions for partial differential
    equations,
    
    .. math::
    
        u = g \hbox{ on } G,
    
    where :math:`u` is the solution to be computed, :math:`g` is a function
    and :math:`G` is a sub domain of the mesh.
    
    A DirichletBC is specified by the function g, the function space
    (trial space) and boundary indicators on (a subset of) the mesh
    boundary.
    
    The boundary indicators may be specified in a number of
    different ways.
    
    The simplest approach is to specify the boundary by a :cpp:class:`SubDomain`
    object, using the inside() function to specify on which facets
    the boundary conditions should be applied. The boundary facets
    will then be searched for and marked *only* on the first call to
    apply. This means that the mesh could be moved after the first
    apply and the boundary markers would still remain intact.
    
    Alternatively, the boundary may be specified by a :cpp:class:`MeshFunction`
    over facets labeling all mesh facets together with a number that
    specifies which facets should be included in the boundary.
    
    The third option is to attach the boundary information to the
    mesh. This is handled automatically when exporting a mesh from
    for example VMTK.
    
    The 'method' variable may be used to specify the type of method
    used to identify degrees of freedom on the boundary. Available
    methods are: topological approach (default), geometric approach,
    and pointwise approach. The topological approach is faster, but
    will only identify degrees of freedom that are located on a
    facet that is entirely on the boundary. In particular, the
    topological approach will not identify degrees of freedom for
    discontinuous elements (which are all internal to the cell). A
    remedy for this is to use the geometric approach. In the
    geometric approach, each dof on each facet that matches the
    boundary condition will be checked. To apply pointwise boundary
    conditions e.g. pointloads, one will have to use the pointwise
    approach. The three possibilities are "topological", "geometric"
    and "pointwise".
    
    Note: when using "pointwise", the boolean argument `on_boundary`
    in SubDomain::inside will always be false.
    
    The 'check_midpoint' variable can be used to decide whether or
    not the midpoint of each facet should be checked when a
    user-defined :cpp:class:`SubDomain` is used to define the domain of the
    boundary condition. By default, midpoints are always checked.
    Note that this variable may be of importance close to corners,
    in which case it is sometimes important to check the midpoint to
    avoid including facets "on the diagonal close" to a corner. This
    variable is also of importance for curved boundaries (like on a
    sphere or cylinder), in which case it is important *not* to
    check the midpoint which will be located in the interior of a
    domain defined relative to a radius.
    
    Note that there may be caching employed in BC computation for
    performance reasons. In particular, applicable DOFs are cached
    by some methods on a first apply(). This means that changing a
    supplied object (defining boundary subdomain) after first use may
    have no effect. But this is implementation and method specific.


    .. cpp:function:: DirichletBC(const FunctionSpace& V, const GenericFunction& g, const SubDomain& sub_domain, std::string method="topological", bool check_midpoint=true)
    
        Create boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain.
            method (std::string)
                Optional argument: A string specifying
                the method to identify dofs.


    .. cpp:function:: DirichletBC(std::shared_ptr<const FunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::shared_ptr<const SubDomain> sub_domain, std::string method="topological", bool check_midpoint=true)
    
        Create boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space
            g (:cpp:class:`GenericFunction`)
                The value
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain
            method (std::string)
                Optional argument: A string specifying
                the method to identify dofs


    .. cpp:function:: DirichletBC(const FunctionSpace& V, const GenericFunction& g, const MeshFunction<std::size_t>& sub_domains, std::size_t sub_domain, std::string method="topological")
    
        Create boundary condition for subdomain specified by index
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domains (:cpp:class:`MeshFunction` <std::size_t>)
                Subdomain markers
            sub_domain (std::size_t)
                The subdomain index (number)
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


    .. cpp:function:: DirichletBC(std::shared_ptr<const FunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::shared_ptr<const MeshFunction<std::size_t> > sub_domains, std::size_t sub_domain, std::string method="topological")
    
        Create boundary condition for subdomain specified by index
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domains (:cpp:class:`MeshFunction` <std::size_t>)
                Subdomain markers
            sub_domain (std::size_t)
                The subdomain index (number)
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


    .. cpp:function:: DirichletBC(const FunctionSpace& V, const GenericFunction& g, std::size_t sub_domain, std::string method="topological")
    
        Create boundary condition for boundary data included in the mesh
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domain (std::size_t)
                The subdomain index (number)
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


    .. cpp:function:: DirichletBC(std::shared_ptr<const FunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::size_t sub_domain, std::string method="topological")
    
        Create boundary condition for boundary data included in the mesh
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domain (std::size_t)
                The subdomain index (number)
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


    .. cpp:function:: DirichletBC(std::shared_ptr<const FunctionSpace> V, std::shared_ptr<const GenericFunction> g, const std::vector<std::size_t>& markers, std::string method="topological")
    
        Create boundary condition for subdomain by boundary markers
        (cells, local facet numbers)
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            markers (std::vector<std::size_t>)
                Subdomain markers (facet index local to process)
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


    .. cpp:function:: DirichletBC(const DirichletBC& bc)
    
        Copy constructor. Either cached DOF data are copied.
        
        *Arguments*
            bc (:cpp:class:`DirichletBC`)
                The object to be copied.


    .. cpp:function:: const DirichletBC& operator= (const DirichletBC& bc)
    
        Assignment operator. Either cached DOF data are assigned.
        
        *Arguments*
            bc (:cpp:class:`DirichletBC`)
                Another DirichletBC object.


    .. cpp:function:: void apply(GenericMatrix& A) const
    
        Apply boundary condition to a matrix
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary condition to.


    .. cpp:function:: void apply(GenericVector& b) const
    
        Apply boundary condition to a vector
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary condition to.


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b) const
    
        Apply boundary condition to a linear system
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary condition to.
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary condition to.


    .. cpp:function:: void apply(GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to vectors for a nonlinear problem
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary conditions to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a linear system for a nonlinear problem
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary conditions to.
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary conditions to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


    .. cpp:function:: void get_boundary_values(Map& boundary_values, std::string method="default") const
    
        Get Dirichlet dofs and values. If a method other than 'pointwise' is
        used in parallel, the map may not be complete for local vertices since
        a vertex can have a bc applied, but the partition might not have a
        facet on the boundary. To ensure all local boundary dofs are marked,
        it is necessary to call gather() on the returned boundary values.
        
        *Arguments*
            boundary_values (std::unordered_map<std::size_t, double>)
                Map from dof to boundary value.
            method (std::string)
                Optional argument: A string specifying which
                method to use.


    .. cpp:function:: void gather(Map& boundary_values) const
    
        Get boundary values from neighbour processes. If a method other than
        "pointwise" is used, this is necessary to ensure all boundary dofs are
        marked on all processes.
        
        *Arguments*
            boundary_values (std::unordered_map<std::size_t, double>)
                Map from dof to boundary value.


    .. cpp:function:: void zero(GenericMatrix& A) const
    
        Make rows of matrix associated with boundary condition zero,
        useful for non-diagonal matrices in a block matrix.
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix


    .. cpp:function:: void zero_columns(GenericMatrix& A, GenericVector& b, double diag_val=0) const
    
        Make columns of matrix associated with boundary condition
        zero, and update a (right-hand side) vector to reflect the
        changes. Useful for non-diagonals.
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix
            b (:cpp:class:`GenericVector`)
                The vector
            diag_val (double)
                This parameter would normally be -1, 0 or 1.


    .. cpp:function:: const std::vector<std::size_t>& markers() const
    
        Return boundary markers
        
        *Returns*
            std::vector<std::pair<std::size_t, std::size_t> >
                Boundary markers (facets stored as pairs of cells and
                local facet numbers).


    .. cpp:function:: std::shared_ptr<const FunctionSpace> function_space() const
    
        Return function space V
        
        *Returns*
            :cpp:class:`FunctionSpace`
                The function space to which boundary conditions are applied.


    .. cpp:function:: std::shared_ptr<const GenericFunction> value() const
    
        Return boundary value g
        
        *Returns*
            :cpp:class:`GenericFunction`
                The boundary values.


    .. cpp:function:: std::shared_ptr<const SubDomain> user_sub_domain() const
    
        Return shared pointer to subdomain
        
        *Returns*
            :cpp:class:`SubDomain`
                Shared pointer to subdomain.


    .. cpp:function:: bool is_compatible(GenericFunction& v) const
    
        Check if given function is compatible with boundary condition
        (checking only vertex values)
        
        *Arguments*
            v (:cpp:class:`GenericFunction`)
                The function to check for compatibility
                with boundary condition.
        
        *Returns*
            bool
                True if compatible.


    .. cpp:function:: void set_value(const GenericFunction& g)
    
        Set value g for boundary condition, domain remains unchanged
        
        *Arguments*
            g (:cpp:class:`GenericFunction`)
                The value.


    .. cpp:function:: void set_value(std::shared_ptr<const GenericFunction> g)
    
        Set value g for boundary condition, domain remains unchanged
        
        *Arguments*
            g (:cpp:class:`GenericFunction`)
                The value.


    .. cpp:function:: void homogenize()
    
        Set value to 0.0


    .. cpp:function:: std::string method() const
    
        Return method used for computing Dirichlet dofs
        
        *Returns*
            std::string
                Method used for computing Dirichlet dofs ("topological",
                "geometric" or "pointwise").


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


.. cpp:class:: LocalData

