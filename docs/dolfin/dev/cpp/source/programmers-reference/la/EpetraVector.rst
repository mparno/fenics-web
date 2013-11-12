
.. Documentation for the header file dolfin/la/EpetraVector.h

.. _programmers_reference_cpp_la_epetravector:

EpetraVector.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EpetraVector

    *Parent class(es)*
    
        * :cpp:class:`GenericVector`
        
    This class provides a simple vector class based on Epetra.
    It is a simple wrapper for an Epetra vector object (Epetra_FEVector)
    implementing the GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the Epetra_FEVector object using the function vec() or vec_ptr()
    and use the standard Epetra interface.


    .. cpp:function:: EpetraVector(std::string type="global")
    
        Create empty vector


    .. cpp:function:: explicit EpetraVector(std::size_t N, std::string type="global")
    
        Create vector of size N


    .. cpp:function:: EpetraVector(const EpetraVector& x)
    
        Copy constructor


    .. cpp:function:: explicit EpetraVector(boost::shared_ptr<Epetra_FEVector> vector)
    
        Create vector view from given Epetra_FEVector pointer


    .. cpp:function:: explicit EpetraVector(const Epetra_BlockMap& map)
    
        Create vector from given Epetra_BlockMap


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<GenericVector> copy() const
    
        Return copy of vector


    .. cpp:function:: void resize(std::size_t N)
    
        Resize vector to size N


    .. cpp:function:: void resize(std::pair<std::size_t, std::size_t> range)
    
        Resize vector with given ownership range


    .. cpp:function:: void resize(std::pair<std::size_t, std::size_t> range, const std::vector<la_index>& ghost_indices)
    
        Resize vector with given ownership range and with ghost values


    .. cpp:function:: bool empty() const
    
        Return true if vector is empty


    .. cpp:function:: std::size_t size() const
    
        Return size of vector


    .. cpp:function:: std::size_t local_size() const
    
        Return size of local vector


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range() const
    
        Return local ownership range of a vector


    .. cpp:function:: bool owns_index(std::size_t i) const
    
        Determine whether global vector index is owned by this process


    .. cpp:function:: void set(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Set block of values


    .. cpp:function:: void add(const double* block, std::size_t m, const dolfin::la_index* rows)
    
        Add block of values


    .. cpp:function:: void get_local(std::vector<double>& values) const
    
        Get all values on local process


    .. cpp:function:: void set_local(const std::vector<double>& values)
    
        Set all values on local process


    .. cpp:function:: void add_local(const Array<double>& values)
    
        Add all values to each entry on local process


    .. cpp:function:: void gather(GenericVector& x, const std::vector<dolfin::la_index>& indices) const
    
        Gather entries into local vector x


    .. cpp:function:: void gather(std::vector<double>& x, const std::vector<dolfin::la_index>& indices) const
    
        Gather entries into x


    .. cpp:function:: void gather_on_zero(std::vector<double>& x) const
    
        Gather all entries into x on process 0


    .. cpp:function:: void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)


    .. cpp:function:: void abs()
    
        Replace all entries in the vector by their absolute values


    .. cpp:function:: double inner(const GenericVector& vector) const
    
        Return inner product with given vector


    .. cpp:function:: double norm(std::string norm_type) const
    
        Return norm of vector


    .. cpp:function:: double min() const
    
        Return minimum value of vector


    .. cpp:function:: double max() const
    
        Return maximum value of vector


    .. cpp:function:: double sum() const
    
        Return sum of values of vector


    .. cpp:function:: double sum(const Array<std::size_t>& rows) const
    
        Return sum of selected rows in vector


    .. cpp:function:: const EpetraVector& operator*= (double a)
    
        Multiply vector by given number


    .. cpp:function:: const EpetraVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise


    .. cpp:function:: const EpetraVector& operator/= (double a)
    
        Divide vector by given number


    .. cpp:function:: const EpetraVector& operator+= (const GenericVector& x)
    
        Add given vector


    .. cpp:function:: const EpetraVector& operator+= (double a)
    
        Add number to all components of a vector


    .. cpp:function:: const EpetraVector& operator-= (const GenericVector& x)
    
        Subtract given vector


    .. cpp:function:: const EpetraVector& operator-= (double a)
    
        Subtract number from all components of a vector


    .. cpp:function:: const EpetraVector& operator= (const GenericVector& x)
    
        Assignment operator


    .. cpp:function:: const EpetraVector& operator= (double a)
    
        Assignment operator


    .. cpp:function:: GenericLinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory


    .. cpp:function:: void reset(const Epetra_BlockMap& map)
    
        Reset Epetra_FEVector


    .. cpp:function:: boost::shared_ptr<Epetra_FEVector> vec() const
    
        Return Epetra_FEVector pointer


    .. cpp:function:: const EpetraVector& operator= (const EpetraVector& x)
    
        Assignment operator


