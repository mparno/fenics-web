
.. Documentation for the header file dolfin/ale/MeshDisplacement.h

.. _programmers_reference_cpp_ale_meshdisplacement:

MeshDisplacement.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshDisplacement

    *Parent class(es)*
    
        * :cpp:class:`Expression`
        
    This class encapsulates the CG1 representation of the
    displacement of a mesh as an Expression. This is particularly
    useful for the displacement returned by mesh smoothers which can
    subsequently be used in evaluating forms. The value rank is 1
    and the value shape is equal to the geometric dimension of the
    mesh.


    .. cpp:function:: MeshDisplacement(const Mesh& mesh)
    
        Create MeshDisplacement of given mesh
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Mesh to be displacement defined on.


    .. cpp:function:: MeshDisplacement(const MeshDisplacement& mesh_displacement)
    
        Copy constructor
        
        *Arguments*
            mesh_displacement (:cpp:class:`MeshDisplacement`)
                Object to be copied.


    .. cpp:function:: Function& operator[] (const std::size_t i)
    
        Extract subfunction
        In python available as MeshDisplacement.sub(i)
        
        *Arguments*
            i (std::size_t)
                Index of subfunction.


    .. cpp:function:: const Function& operator[] (const std::size_t i) const
    
        Extract subfunction. Const version
        
        *Arguments*
            i (std::size_t)
                Index of subfunction.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Evaluate at given point in given cell.
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values at the point.
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.
            cell (ufc::cell)
                The cell which contains the given point.


    .. cpp:function:: void compute_vertex_values(std::vector<double>& vertex_values, const Mesh& mesh) const
    
        Compute values at all mesh vertices.
        
        *Arguments*
            vertex_values (:cpp:class:`Array` <double>)
                The values at all vertices.
            mesh (:cpp:class:`Mesh`)
                The mesh.


