
.. Documentation for the header file dolfin/fem/Form.h

.. _programmers_reference_cpp_fem_form:

Form.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Form

    *Parent class(es)*
    
        * :cpp:class:`Hierarchical<Form>`
        
    A note on the order of trial and test spaces: FEniCS numbers
    argument spaces starting with the leading dimension of the
    corresponding tensor (matrix). In other words, the test space is
    numbered 0 and the trial space is numbered 1. However, in order
    to have a notation that agrees with most existing finite element
    literature, in particular
    
        a = a(u, v)
    
    the spaces are numbered from right to left:
    
       a : V_1 x V_0 --> R
    
    This is reflected in the ordering of the spaces that should be
    supplied to generated subclasses. In particular, when a bilinear
    form is initialized, it should be initialized as
    
       a(V1, V0)
    
    where V1 is the trial space and V0 is the test space. However,
    when a form is initialized by a list of argument spaces (the
    variable 'function_spaces' in the constructors below, the list
    of spaces should start with space number 0 (the test space) and
    then space number 1 (the trial space).


    .. cpp:function:: Form(dolfin::uint rank, dolfin::uint num_coefficients)
    
        Create form of given rank with given number of coefficients


    .. cpp:function:: Form(boost::shared_ptr<const ufc::form> ufc_form, std::vector<boost::shared_ptr<const FunctionSpace> > function_spaces, std::vector<boost::shared_ptr<const GenericFunction> > coefficients)
    
        Create form (shared data)


    .. cpp:function:: uint rank() const
    
        Return rank of form (bilinear form = 2, linear form = 1,
        functional = 0, etc)


    .. cpp:function:: uint num_coefficients() const
    
        Return number of coefficients


    .. cpp:function:: std::vector<uint> coloring(uint entity_dim) const
    
        Return coloring type for colored (multi-threaded) assembly of form
        over a mesh entity of a given dimension


    .. cpp:function:: void set_mesh(boost::shared_ptr<const Mesh> mesh)
    
        Set mesh, necessary for functionals when there are no function spaces


    .. cpp:function:: const Mesh& mesh() const
    
        Extract common mesh from form


    .. cpp:function:: boost::shared_ptr<const Mesh> mesh_shared_ptr() const
    
        Return mesh shared pointer (if any)


    .. cpp:function:: boost::shared_ptr<const FunctionSpace> function_space(uint i) const
    
        Return function space for given argument


    .. cpp:function:: std::vector<boost::shared_ptr<const FunctionSpace> > function_spaces() const
    
        Return function spaces for arguments


    .. cpp:function:: void set_coefficient(uint i, boost::shared_ptr<const GenericFunction> coefficient)
    
        Set coefficient with given number (shared pointer version)


    .. cpp:function:: void set_coefficient(std::string name, boost::shared_ptr<const GenericFunction> coefficient)
    
        Set coefficient with given name (shared pointer version)


    .. cpp:function:: void set_coefficients(std::map<std::string, boost::shared_ptr<const GenericFunction> > coefficients)
    
        Set all coefficients in given map, possibly a subset
        (shared pointer version)


    .. cpp:function:: boost::shared_ptr<const GenericFunction> coefficient(uint i) const
    
        Return coefficient with given number


    .. cpp:function:: boost::shared_ptr<const GenericFunction> coefficient(std::string name) const
    
        Return coefficient with given name


    .. cpp:function:: std::vector<boost::shared_ptr<const GenericFunction> > coefficients() const
    
        Return all coefficients


    .. cpp:function:: dolfin::uint coefficient_number(const std::string & name) const
    
        Return the number of the coefficient with this name


    .. cpp:function:: std::string coefficient_name(dolfin::uint i) const
    
        Return the name of the coefficient with this number


    .. cpp:function:: boost::shared_ptr<const MeshFunction<uint> > cell_domains_shared_ptr() const
    
        Return cell domains (pointer may be zero if no domains have been specified)


    .. cpp:function:: boost::shared_ptr<const MeshFunction<uint> > exterior_facet_domains_shared_ptr() const
    
        Return exterior facet domains (pointer may be zero if no domains have been specified)


    .. cpp:function:: boost::shared_ptr<const MeshFunction<uint> > interior_facet_domains_shared_ptr() const
    
        Return interior facet domains (pointer may be zero if no domains have been specified)


    .. cpp:function:: void set_cell_domains(boost::shared_ptr<const MeshFunction<unsigned int> > cell_domains)
    
        Set cell domains


    .. cpp:function:: void set_exterior_facet_domains(boost::shared_ptr<const MeshFunction<unsigned int> > exterior_facet_domains)
    
        Set exterior facet domains


    .. cpp:function:: void set_interior_facet_domains(boost::shared_ptr<const MeshFunction<unsigned int> > interior_facet_domains)
    
        Set interior facet domains


    .. cpp:function:: const ufc::form& ufc_form() const
    
        Return UFC form


    .. cpp:function:: boost::shared_ptr<const ufc::form> ufc_form_shared_ptr() const
    
        Return UFC form shared pointer


    .. cpp:function:: void check() const
    
        Check function spaces and coefficients


