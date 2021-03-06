
.. Documentation for the header file dolfin/intersection/IntersectionOperator.h

.. _programmers_reference_cpp_intersection_intersectionoperator:

IntersectionOperator.h
======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: IntersectionOperator

    .. cpp:function:: IntersectionOperator(const Mesh& _mesh, const std::string& kernel_type = "SimpleCartesian")
    
        Create intersection detector for the mesh \em mesh.
        @param kernel_type The CGAL geometric kernel is used to compute predicates,
        intersections and such. Depending on this choice the kernel
        (kernel_type = "ExcactPredicates") can compute predicates excactly
        (without roundoff error) or only approximately (default, kernel_type =
        "SimpleCartesian").


    .. cpp:function:: void all_intersected_entities(const Point & point, std::set<uint>& ids_result) const
    
        Compute all id of all cells which are intersects by a \em point.
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.


    .. cpp:function:: void all_intersected_entities(const std::vector<Point>& points, std::set<uint>& ids_result) const
    
        Compute all id of all cells which are intersects any point in \em points.
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.


    .. cpp:function:: void all_intersected_entities(const MeshEntity & entity, std::vector<uint> & ids_result) const
    
        Compute all id of all cells which are intersects by a \em entity.
        \param[out] ids_result The ids of the intersected entities are saved in a vector.
        This allows is more efficent than using a set and allows a map between
        the (external) cell and the intersected cell of the mesh. If you
        are only interested in intersection with a list of cells without caring about which
        cell what intersected by which one, use
        void IntersectionOperator::all_intersected_entities(const std::vector<Cell> &, std::set<uint> &) const;
        @internal
        @todo This function has to improved: 1) it requires the object the
        mesh is to be cut with to be another mesh entitiy instead of being just a
        kind of geometric object. 2) Requires a runtime switch 3) would require a
        implementation for each geometric  primitive if they have no common base
        class.


    .. cpp:function:: void all_intersected_entities(const std::vector<MeshEntity> & entities, std::set<uint> & ids_result) const
    
        Compute all id of all cells which are intersects by any of the entities in \em entities. This
        \param[out] ids_result The ids of the intersected set are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.


    .. cpp:function:: void all_intersected_entities(const Mesh& another_mesh, std::set<uint>& ids_result) const
    
        Compute all id of all cells which are intersects by the given mesh \em another_mesh;
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.


    .. cpp:function:: int any_intersected_entity(const Point& point) const
    
        Computes only the first id of the entity, which contains the point. Returns -1 if no cell is intersected.
        @internal @remark This makes the function evaluation significantly faster.


    .. cpp:function:: Point closest_point(const Point & point) const
    
        Computes the point inside the mesh which is closest to the point query.


    .. cpp:function:: dolfin::uint closest_cell(const Point & point) const
    
        Computes the index of the cell inside the mesh which are closest to the point query.


    .. cpp:function:: std::pair<Point,uint> closest_point_and_cell(const Point & point) const
    
        Computes the point inside the mesh and the corresponding cell index
        that are closest to the point query.


    .. cpp:function:: void reset_kernel(const std::string& kernel_type  = "SimpleCartesian")
    
        Rebuilds the underlying search structure from scratch and uses
        the kernel kernel_type underlying CGAL Geometry kernel.


    .. cpp:function:: void clear()
    
        Clears search structure. Should be used if the mesh has changed


