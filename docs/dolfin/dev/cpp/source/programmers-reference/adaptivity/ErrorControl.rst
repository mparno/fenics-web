
.. Documentation for the header file dolfin/adaptivity/ErrorControl.h

.. _programmers_reference_cpp_adaptivity_errorcontrol:

ErrorControl.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ErrorControl

    *Parent class(es)*
    
        * :cpp:class:`Hierarchical<ErrorControl>`
        
        * :cpp:class:`Variable`
        
    (Goal-oriented) Error Control class.
    The notation used here follows the notation in "Automated
    goal-oriented error control I: stationary variational problems",
    ME Rognes and A Logg, 2010-2011.


    .. cpp:function:: ErrorControl(boost::shared_ptr<Form> a_star, boost::shared_ptr<Form> L_star, boost::shared_ptr<Form> residual, boost::shared_ptr<Form> a_R_T, boost::shared_ptr<Form> L_R_T, boost::shared_ptr<Form> a_R_dT, boost::shared_ptr<Form> L_R_dT, boost::shared_ptr<Form> eta_T, bool is_linear)
    
        Create error control object
        
        *Arguments*
            a_star (:cpp:class:`Form`)
               the bilinear form for the dual problem
            L_star (:cpp:class:`Form`)
               the linear form for the dual problem
            residual (:cpp:class:`Form`)
               a functional for the residual (error estimate)
            a_R_T (:cpp:class:`Form`)
               the bilinear form for the strong cell residual problem
            L_R_T (:cpp:class:`Form`)
               the linear form for the strong cell residual problem
            a_R_dT (:cpp:class:`Form`)
               the bilinear form for the strong facet residual problem
            L_R_dT (:cpp:class:`Form`)
               the linear form for the strong facet residual problem
            eta_T (:cpp:class:`Form`)
               a linear form over DG_0 for error indicators
            is_linear (bool)
               true iff primal problem is linear


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values:


    .. cpp:function:: double estimate_error(const Function& u, const std::vector<boost::shared_ptr<const DirichletBC> > bcs)
    
        Estimate the error relative to the goal M of the discrete
        approximation 'u' relative to the variational formulation by
        evaluating the weak residual at an approximation to the dual
        solution.
        
        *Arguments*
            u (:cpp:class:`Function`)
               the primal approximation
        
            bcs (std::vector<:cpp:class:`DirichletBC`>)
                the primal boundary conditions
        
        *Returns*
            double
                error estimate


    .. cpp:function:: void compute_indicators(MeshFunction<double>& indicators, const Function& u)
    
        Compute error indicators
        
        *Arguments*
            indicators (:cpp:class:`MeshFunction` <double>)
                the error indicators (to be computed)
        
            u (:cpp:class:`Function`)
                the primal approximation


    .. cpp:function:: void residual_representation(Function& R_T, SpecialFacetFunction& R_dT, const Function& u)
    
        Compute strong representation (strong cell and facet
        residuals) of the weak residual.
        
        *Arguments*
            R_T (:cpp:class:`Function`)
                the strong cell residual (to be computed)
        
            R_dT (:cpp:class:`SpecialFacetFunction`)
                the strong facet residual (to be computed)
        
            u (:cpp:class:`Function`)
                the primal approximation


    .. cpp:function:: void compute_cell_residual(Function& R_T, const Function& u)
    
        Compute representation for the strong cell residual
        from the weak residual
        
        *Arguments*
            R_T (:cpp:class:`Function`)
                the strong cell residual (to be computed)
        
            u (:cpp:class:`Function`)
                the primal approximation


    .. cpp:function:: void compute_facet_residual(SpecialFacetFunction& R_dT, const Function& u, const Function& R_T)
    
        Compute representation for the strong facet residual from the
        weak residual and the strong cell residual
        
        *Arguments*
            R_dT (:cpp:class:`SpecialFacetFunction`)
                the strong facet residual (to be computed)
        
            u (:cpp:class:`Function`)
                the primal approximation
        
            R_T (:cpp:class:`Function`)
                the strong cell residual


    .. cpp:function:: void compute_dual(Function& z, const std::vector<boost::shared_ptr<const DirichletBC> > bcs)
    
        Compute dual approximation defined by dual variational
        problem and dual boundary conditions given by homogenized primal
        boundary conditions.
        
        *Arguments*
            z (:cpp:class:`Function`)
                the dual approximation (to be computed)
        
            bcs (std::vector<:cpp:class:`DirichletBC`>)
                the primal boundary conditions


    .. cpp:function:: void compute_extrapolation(const Function& z, const std::vector<boost::shared_ptr<const DirichletBC> > bcs)
    
        Compute extrapolation with boundary conditions
        
        *Arguments*
            z (:cpp:class:`Function`)
                the extrapolated function (to be computed)
        
            bcs (std::vector<:cpp:class:`DirichletBC`>)
                the dual boundary conditions


