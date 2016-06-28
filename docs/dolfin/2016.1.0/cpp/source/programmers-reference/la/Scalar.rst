
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


    .. cpp:function:: Scalar(MPI_Comm comm)
    
        Create zero scalar


    .. cpp:function:: void init(const TensorLayout& tensor_layout)
    
        Initialize zero tensor using sparsity pattern


    .. cpp:function:: bool empty() const
    
        Return true if empty


    .. cpp:function:: std::size_t rank() const
    
        Return tensor rank (number of dimensions)


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: std::pair<std::int64_t, std::int64_t> local_range(std::size_t dim) const
    
        Return local ownership range


    .. cpp:function:: void get(double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows) const
    
        Get block of values


    .. cpp:function:: void set(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Set block of values using global indices


    .. cpp:function:: void set_local(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Set block of values using local indices


    .. cpp:function:: void add(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Add block of values using global indices


    .. cpp:function:: void add_local(const double* block, const dolfin::la_index* num_rows, const dolfin::la_index * const * rows)
    
        Add block of values using local indices


    .. cpp:function:: void add(const double* block, const std::vector<ArrayView<const dolfin::la_index>>& rows)
    
        Add block of values using global indices


    .. cpp:function:: void add_local(const double* block, const std::vector<ArrayView<const dolfin::la_index>>& rows)
    
        Add block of values using local indices


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::shared_ptr<Scalar> copy() const
    
        Return copy of scalar


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return a factory for the default linear algebra backend


    .. cpp:function:: double get_scalar_value() const
    
        Get final value (assumes prior apply(), not part of
        GenericTensor interface)


    .. cpp:function:: void add_local_value(double value)
    
        Add to local increment (added for testing, remove if we add a
        better way from python)


