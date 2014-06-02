
.. Documentation for the header file dolfin/multistage/MultiStageScheme.h

.. _programmers_reference_cpp_multistage_multistagescheme:

MultiStageScheme.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiStageScheme

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: MultiStageScheme(std::vector<std::vector<std::shared_ptr<const Form> > > stage_forms, std::shared_ptr<const Form> last_stage, std::vector<std::shared_ptr<Function> > stage_solutions, std::shared_ptr<Function> u, std::shared_ptr<Constant> t, std::shared_ptr<Constant> dt, std::vector<double> dt_stage_offset, std::vector<int> jacobian_indices, unsigned int order, const std::string name, const std::string human_form)
    
        Constructor
        FIXME: This constructor is a MESS. Needs clean up...


    .. cpp:function:: MultiStageScheme(std::vector<std::vector<std::shared_ptr<const Form> > > stage_forms, std::shared_ptr<const Form> last_stage, std::vector<std::shared_ptr<Function> > stage_solutions, std::shared_ptr<Function> u, std::shared_ptr<Constant> t, std::shared_ptr<Constant> dt, std::vector<double> dt_stage_offset, std::vector<int> jacobian_indices, unsigned int order, const std::string name, const std::string human_form, std::vector<const DirichletBC* > bcs)
    
        Constructor with Boundary conditions


    .. cpp:function:: std::vector<std::vector<std::shared_ptr<const Form> > >& stage_forms()
    
        Return the stages


    .. cpp:function:: std::shared_ptr<const Form> last_stage()
    
        Return the last stage


    .. cpp:function:: std::vector<std::shared_ptr<Function> >& stage_solutions()
    
        Return stage solutions


    .. cpp:function:: std::shared_ptr<Function> solution()
    
        Return solution variable


    .. cpp:function:: std::shared_ptr<const Function> solution() const
    
        Return solution variable (const version)


    .. cpp:function:: std::shared_ptr<Constant> t()
    
        Return local time


    .. cpp:function:: std::shared_ptr<Constant> dt()
    
        Return local timestep


    .. cpp:function:: const std::vector<double>& dt_stage_offset() const
    
        Return local timestep


    .. cpp:function:: unsigned int order() const
    
        Return the order of the scheme


    .. cpp:function:: std::vector<const DirichletBC* > bcs() const
    
        Return boundary conditions


    .. cpp:function:: bool implicit(unsigned int stage) const
    
        Return true if stage is implicit


    .. cpp:function:: bool implicit() const
    
        Return true if the whole scheme is implicit


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


