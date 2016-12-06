
.. Documentation for the header file dolfin/mesh/MeshQuality.h

.. _programmers_reference_cpp_mesh_meshquality:

MeshQuality.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshQuality

    The class provides functions to quantify mesh quality


    .. cpp:function:: static CellFunction<double> radius_ratios(std::shared_ptr<const Mesh> mesh)
    
        Compute the radius ratio for all cells.
        
        *Returns*
            CellFunction<double>
                The cell radius ratio radius ratio geometric_dimension *
                * inradius / circumradius (geometric_dimension
                is normalization factor). It has range zero to one.
                Zero indicates a degenerate element.
        
        *Example*
            .. note::
        
                std::shared_ptr<Mesh> mesh(new UnitCubeMesh(4, 4, 4));
                CellFunction<double> = MeshQuality::radius_ratio(mesh);


    .. cpp:function:: static std::pair<double, double> radius_ratio_min_max(const Mesh& mesh)
    
        Compute the minimum and maximum radius ratio of cells
        (across all processes)
        
        *Returns*
            std::pair<double, double>
                The [minimum, maximum] cell radii ratio (geometric_dimension *
                * inradius / circumradius, geometric_dimension
                is normalization factor). It has range zero to one.
                Zero indicates a degenerate element.
        
        *Example*
            .. note::
        
                Mesh  UnitCubeMesh(4, 4, 4);
                std::pair<double, double> ratios
                   = MeshQuality::radius_ratio_min_max(mesh);
                double min_ratio = ratios.first;
                double max_ratio = ratios.second;


    .. cpp:function:: static std::pair<std::vector<double>, std::vector<double> > radius_ratio_histogram_data(const Mesh& mesh, std::size_t num_bins = 50)
    
        Create (ratio, number of cells) data for creating a histogram
        of cell quality


    .. cpp:function:: static std::string radius_ratio_matplotlib_histogram(const Mesh& mesh, std::size_t num_intervals = 50)
    
        Create Matplotlib string to plot cell quality histogram


    .. cpp:function:: static void dihedral_angles(const Cell& cell, std::vector<double>& dh_angle)
    
        Get internal dihedral angles of a tetrahedral cell


    .. cpp:function:: static std::pair<double, double> dihedral_angles_min_max(const Mesh& mesh)
    
        Get internal minimum and maximum dihedral angles of a 3D mesh


    .. cpp:function:: static std::pair<std::vector<double>, std::vector<double> > dihedral_angles_histogram_data(const Mesh& mesh, std::size_t num_bins = 100)
    
        Create (dihedral angles, number of cells) data for creating a histogram
        of dihedral


    .. cpp:function:: static std::string dihedral_angles_matplotlib_histogram(const Mesh& mesh, std::size_t num_intervals = 100)
    
        Create Matplotlib string to plot dihedral angles quality histogram


