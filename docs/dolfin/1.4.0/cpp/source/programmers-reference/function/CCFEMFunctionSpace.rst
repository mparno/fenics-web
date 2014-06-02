
.. Documentation for the header file dolfin/function/CCFEMFunctionSpace.h

.. _programmers_reference_cpp_function_ccfemfunctionspace:

CCFEMFunctionSpace.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMFunctionSpace

    This class represents a cut and composite finite element
    function space (CCFEM) defined on one or more possibly
    intersecting meshes.
    
    A CCFEM function space may be created from a set of standard
    function spaces by repeatedly calling add(), followed by a call
    to build(). Note that a CCFEM function space is not useful and
    its data structures are empty until build() has been called.


    .. cpp:function:: CCFEMFunctionSpace()
    
        Create empty CCFEM function space


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of the CCFEM function space
        
        *Returns*
            std::size_t
                The dimension of the CCFEM function space.


    .. cpp:function:: std::shared_ptr<const CCFEMDofMap> dofmap() const
    
        Return CCFEM dofmap
        
        *Returns*
            :cpp:class:`CCFEMDofMap`
                The dofmap.


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number function spaces (parts) of the CCFEM function space
        
        *Returns*
            std::size_t
                The number of function spaces (parts) of the CCFEM function space.


    .. cpp:function:: std::shared_ptr<const FunctionSpace> part(std::size_t i) const
    
        Return function space (part) number i
        
        *Arguments*
            i (std::size_t)
                The part number
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Function space (part) number i


    .. cpp:function:: const std::vector<unsigned int>& uncut_cells(std::size_t part) const
    
        Return the list of uncut cells for given part. The uncut cells
        are defined as all cells that don't collide with any cells in
        any other part with higher part number.
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of uncut cell indices for given part


    .. cpp:function:: const std::vector<unsigned int>& cut_cells(std::size_t part) const
    
        Return the list of cut cells for given part. The cut cells are
        defined as all cells that collide with the boundary of any
        part with higher part number.
        
        FIXME: Figure out whether this makes sense; a cell may collide
        with the boundary of part j but may still be covered
        completely by the domain of part j + 1. Possible solution is
        to for each part i check overlapping parts starting from the
        top and working back down to i + 1.
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of cut cell indices for given part


    .. cpp:function:: const std::vector<unsigned int>& covered_cells(std::size_t part) const
    
        Return the list of covered cells for given part. The covered
        cells are defined as all cells that collide with the domain of
        any part with higher part number, but not with the boundary of
        that part; in other words cells that are completely covered by
        any other part (and which therefore are inactive).
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of covered cell indices for given part


    .. cpp:function:: const std::map<unsigned int, std::vector<std::pair<std::size_t, unsigned int> > >& collision_map_cut_cells(std::size_t part) const
    
        Return the collision map for cut cells of the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::vector<std::pair<std::size_t, unsigned int> > >
                A map from cell indices of cut cells to a list of
                cutting cells. Each cutting cell is represented as a
                pair (part_number, cutting_cell_index).


    .. cpp:function:: const std::map<unsigned int, std::pair<std::vector<double>, std::vector<double> > > & quadrature_rule_cut_cells(std::size_t part) const
    
        Return quadrature rules for cut cells of the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::pair<std::vector<double>, std::vector<double> > >
                A map from cell indices of cut cells to a quadrature
                rules. Each quadrature rule is represented as a pair
                of an array of quadrature weights and a corresponding
                flattened array of quadrature points.


    .. cpp:function:: void add(std::shared_ptr<const FunctionSpace> function_space)
    
        Add function space (shared pointer version)
        
        *Arguments*
            function_space (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: void add(const FunctionSpace& function_space)
    
        Add function space (reference version)
        
        *Arguments*
            function_space (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: void build()
    
        Build CCFEM function space


