
.. Documentation for the header file dolfin/la/GenericTensor.h

.. _programmers_reference_cpp_la_generictensor:

GenericTensor.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericTensor

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraObject`
        
    This class defines a common interface for arbitrary rank tensors.


    .. cpp:function:: void init(const TensorLayout& tensor_layout) = 0
    
        Initialize zero tensor using tensor layout


    .. cpp:function:: std::size_t rank() const = 0
    
        Return tensor rank (number of dimensions)


    .. cpp:function:: std::size_t size(std::size_t dim) const = 0
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const = 0
    
        Return local ownership range


    .. cpp:function:: void get(double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows) const = 0
    
        Get block of values


    .. cpp:function:: void set(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows) = 0
    
        Set block of values


    .. cpp:function:: void add(const double* block, const std::vector<const std::vector<dolfin::la_index>* >& rows) = 0
    
        Add block of values


    .. cpp:function:: void add(const double* block, const std::vector<std::vector<dolfin::la_index> >& rows) = 0
    
        Add block of values


    .. cpp:function:: void add(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows) = 0
    
        Add block of values


    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const = 0
    
        Return linear algebra backend factory


