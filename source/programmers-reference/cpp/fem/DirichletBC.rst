.. Documentation for the header file dolfin/fem/DirichletBC.h

.. _programmers_reference_cpp_fem_dirichletbc:

DirichletBC.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: DirichletBC

    *Parent class*
    
        * :cpp:class:`BoundaryCondition`
        
    This class specifies the interface for setting (strong)
    Dirichlet boundary conditions for partial differential
    equations,
    
       u = g on G,
    
    where u is the solution to be computed, g is a function
    and G is a sub domain of the mesh.
    
    A DirichletBC is specified by the function g, the function space
    (trial space) and boundary indicators on (a subset of) the mesh
    boundary.
    
    The boundary indicators may be specified in a number of
    different ways.
    
    The simplest approach is to specify the boundary by a SubDomain
    object, using the inside() function to specify on which facets
    the boundary conditions should be applied.
    
    Alternatively, the boundary may be specified by a MeshFunction
    labeling all mesh facets together with a number that specifies
    which facets should be included in the boundary.
    
    The third option is to attach the boundary information to the
    mesh. This is handled automatically when exporting a mesh from
    for example VMTK.
    
    The BCMethod variable may be used to specify the type of method
    used to identify degrees of freedom on the boundary. Available
    methods are: topological approach (default), geometric approach,
    and pointwise approach. The topological approach is faster, but
    will only identify degrees of freedom that are located on a
    facet that is entirely on the boundary. In particular, the
    topological approach will not identify degrees of freedom for
    discontinuous elements (which are all internal to the cell).  A
    remedy for this is to use the geometric approach. To apply
    pointwise boundary conditions e.g. pointloads, one will have to
    use the pointwise approach which in turn is the slowest of the
    three possible methods.  The three possibilties are
    "topological", "geometric" and "pointwise".
    This class specifies the interface for setting (strong)

    .. cpp:function:: DirichletBC(boost::shared_ptr<const FunctionSpace> V,
                                  boost::shared_ptr<const GenericFunction> g,
                                  boost::shared_ptr<const SubDomain> sub_domain,
                                  std::string method="topological")
    
        Create boundary condition for subdomain

    .. cpp:function:: DirichletBC(boost::shared_ptr<const FunctionSpace> V,
                                  boost::shared_ptr<const GenericFunction> g,
                                  const MeshFunction<uint>& sub_domains, uint sub_domain,
                                  std::string method="topological")
    
        Create boundary condition for subdomain specified by index

    .. cpp:function:: DirichletBC(boost::shared_ptr<const FunctionSpace> V,
                                  boost::shared_ptr<const GenericFunction> g,
                                  uint sub_domain,
                                  std::string method="topological")
    
        Create boundary condition for boundary data included in the mesh

    .. cpp:function:: DirichletBC(const DirichletBC& bc)
    
        Copy constructor

    .. cpp:function:: DirichletBC(const FunctionSpace& V,
                                  const GenericFunction& g,
                                  const MeshFunction<uint>& sub_domains, uint sub_domain,
                                  std::string method="topological")
    
        Create boundary condition for subdomain specified by index

    .. cpp:function:: DirichletBC(const FunctionSpace& V,
                                  const GenericFunction& g,
                                  uint sub_domain,
                                  std::string method="topological")
    
        Create boundary condition for boundary data included in the mesh

    .. cpp:function:: DirichletBC(const FunctionSpace& V,
                       const GenericFunction& g,
                       const SubDomain& sub_domain,
                       std::string method="topological")
    
        Create boundary condition for subdomain

    .. cpp:function:: bool is_compatible(GenericFunction& v) const
    
        Check if given function is compatible with boundary condition (checking only vertex values)

    .. cpp:function:: boost::shared_ptr<const GenericFunction> value_ptr()
    
        Return shared pointer to boundary value g
        Testing multiline comment

    .. cpp:function:: const DirichletBC& operator= (const DirichletBC& bc)
    
        Assignment operator

    .. cpp:function:: const GenericFunction& value()
    
        Return boundary value g

    .. cpp:function:: const std::vector<std::pair<uint, uint> >& markers()
    
        Return boundary markers (facets stored as pairs of cells and local facet numbers)

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: void apply(GenericMatrix& A) const
    
        Apply boundary condition to a matrix

    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b) const
    
        Apply boundary condition to a linear system

    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a linear system for a nonlinear problem

    .. cpp:function:: void apply(GenericVector& b) const
    
        Apply boundary condition to a vector

    .. cpp:function:: void apply(GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a vector for a nonlinear problem

    .. cpp:function:: void get_bc(uint* indicators, double* values) const
    
        Get Dirichlet values and indicators

    .. cpp:function:: void set_value(boost::shared_ptr<const GenericFunction> g)
    
        Set value g for boundary condition, domain remains unchanged

    .. cpp:function:: void set_value(const GenericFunction& g)
    
        Set value g for boundary condition, domain remains unchanged

    .. cpp:function:: void zero(GenericMatrix& A) const
    
        Make row associated with boundary conditions zero, useful for non-diagonal matrices in a block matrix.

    .. cpp:function:: ~DirichletBC()
    
        Destructor

