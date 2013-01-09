
.. Documentation for the header file dolfin/la/Scalar.h

.. _programmers_reference_cpp_la_scalar:

Scalar.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Scalar

    *Parent class(es)*
    
        * :cpp:class:`GenericTensor`
        
    This class represents a real-valued scalar quantity and
    implements the GenericTensor interface for scalars.


    .. cpp:function:: Scalar()
    
        Create zero scalar


    .. cpp:function:: void resize(std::size_t rank, const std::size_t* dims)
    
        Resize tensor to given dimensions


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using sparsity pattern


    .. cpp:function:: std::size_t rank() const
    
        Return tensor rank (number of dimensions)


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: void get(double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows) const
    
        Get block of values


    .. cpp:function:: void set(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Set block of values


    .. cpp:function:: void add(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Add block of values


    .. cpp:function:: void add(const double* block, const std::vector<const std::vector<dolfin::la_index>* >& rows)
    
        Add block of values


    .. cpp:function:: void add(const double* block, const std::vector<std::vector<dolfin::la_index> >& rows)
    
        Add block of values


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<Scalar> copy() const
    
        Return copy of scalar


    .. cpp:function:: operator double() const
    
        Cast to double


    .. cpp:function:: const Scalar& operator= (double value)
    
        Assignment from double


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return a factory for the default linear algebra backend


    .. cpp:function:: double getval() const
    
        Get value


