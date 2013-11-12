
.. Documentation for the header file dolfin/la/GenericVector.h

.. _programmers_reference_cpp_la_genericvector:

GenericVector.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericVector

    *Parent class(es)*
    
        * :cpp:class:`GenericTensor`
        
    This class defines a common interface for vectors.


    .. cpp:function:: void resize(std::size_t rank, const std::size_t* dims)
    
        Resize tensor with given dimensions


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


    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<GenericVector> copy() const = 0
    
        Return copy of vector


    .. cpp:function:: void resize(std::size_t N) = 0
    
        Resize vector to global size N


    .. cpp:function:: void resize(std::pair<std::size_t, std::size_t> range) = 0
    
        Resize vector with given ownership range


    .. cpp:function:: void resize(std::pair<std::size_t, std::size_t> range, const std::vector<la_index>& ghost_indices) = 0
    
        Resize vector with given ownership range and with ghost values


    .. cpp:function:: bool empty() const = 0
    
        Return true if empty


    .. cpp:function:: std::size_t size() const = 0
    
        Return global size of vector


    .. cpp:function:: std::size_t local_size() const = 0
    
        Return local size of vector


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range() const = 0
    
        Return local ownership range of a vector


    .. cpp:function:: bool owns_index(std::size_t i) const = 0
    
        Determine whether global vector index is owned by this process


    .. cpp:function:: void get(double* block, std::size_t m, const dolfin::la_index* rows) const
    
        Get block of values (values may live on any process)


    .. cpp:function:: void get_local(double* block, std::size_t m, const dolfin::la_index* rows) const = 0
    
        Get block of values (values must all live on the local process)


    .. cpp:function:: void set(const double* block, std::size_t m, const dolfin::la_index* rows) = 0
    
        Set block of values


    .. cpp:function:: void add(const double* block, std::size_t m, const dolfin::la_index* rows) = 0
    
        Add block of values


    .. cpp:function:: void get_local(std::vector<double>& values) const = 0
    
        Get all values on local process


    .. cpp:function:: void set_local(const std::vector<double>& values) = 0
    
        Set all values on local process


    .. cpp:function:: void add_local(const Array<double>& values) = 0
    
        Add values to each entry on local process


    .. cpp:function:: void gather(GenericVector& x, const std::vector<dolfin::la_index>& indices) const = 0
    
        Gather entries into local vector x


    .. cpp:function:: void gather(std::vector<double>& x, const std::vector<dolfin::la_index>& indices) const = 0
    
        Gather entries into x


    .. cpp:function:: void gather_on_zero(std::vector<double>& x) const = 0
    
        Gather all entries into x on process 0


    .. cpp:function:: void axpy(double a, const GenericVector& x) = 0
    
        Add multiple of given vector (AXPY operation)


    .. cpp:function:: void abs() = 0
    
        Replace all entries in the vector by their absolute values


    .. cpp:function:: double inner(const GenericVector& x) const = 0
    
        Return inner product with given vector


    .. cpp:function:: double norm(std::string norm_type) const = 0
    
        Return norm of vector


    .. cpp:function:: double min() const = 0
    
        Return minimum value of vector


    .. cpp:function:: double max() const = 0
    
        Return maximum value of vector


    .. cpp:function:: double sum() const = 0
    
        Return sum of vector


    .. cpp:function:: double sum(const Array<std::size_t>& rows) const = 0
    
        Return sum of selected rows in vector. Repeated entries are
        only summed once.


    .. cpp:function:: const GenericVector& operator*= (double a) = 0
    
        Multiply vector by given number


    .. cpp:function:: const GenericVector& operator*= (const GenericVector& x) = 0
    
        Multiply vector by another vector pointwise


    .. cpp:function:: const GenericVector& operator/= (double a) = 0
    
        Divide vector by given number


    .. cpp:function:: const GenericVector& operator+= (const GenericVector& x) = 0
    
        Add given vector


    .. cpp:function:: const GenericVector& operator+= (double a) = 0
    
        Add number to all components of a vector


    .. cpp:function:: const GenericVector& operator-= (const GenericVector& x) = 0
    
        Subtract given vector


    .. cpp:function:: const GenericVector& operator-= (double a) = 0
    
        Subtract number from all components of a vector


    .. cpp:function:: const GenericVector& operator= (const GenericVector& x) = 0
    
        Assignment operator


    .. cpp:function:: const GenericVector& operator= (double a) = 0
    
        Assignment operator


    .. cpp:function:: const double* data() const
    
        Return pointer to underlying data (const version)


    .. cpp:function:: double* data()
    
        Return pointer to underlying data


    .. cpp:function:: void update_ghost_values()
    
        Update ghost values


    .. cpp:function:: double operator[] (dolfin::la_index i) const
    
        Get value of given entry


    .. cpp:function:: double getitem(dolfin::la_index i) const
    
        Get value of given entry


    .. cpp:function:: void setitem(dolfin::la_index i, double value)
    
        Set given entry to value. apply("insert") should be called
        before using using the object.


