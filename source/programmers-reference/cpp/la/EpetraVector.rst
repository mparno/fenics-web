.. Documentation for the header file dolfin/la/EpetraVector.h

.. _programmers_reference_cpp_la_epetravector:

EpetraVector.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: EpetraVector

    *Parent class*
    
        * :cpp:class:`GenericVector`
        
    This class provides a simple vector class based on Epetra.
    It is a simple wrapper for an Epetra vector object (Epetra_FEVector)
    implementing the GenericVector interface.
    
    The interface is intentionally simple. For advanced usage,
    access the Epetra_FEVector object using the function vec() or vec_ptr()
    and use the standard Epetra interface.

    .. cpp:function:: EpetraVector(const EpetraVector& x)
    
        Copy constructor

    .. cpp:function:: EpetraVector(std::string type="global")
    
        Create empty vector

    .. cpp:function:: boost::shared_ptr<Epetra_FEVector> vec() const
    
        Return Epetra_FEVector pointer

    .. cpp:function:: const EpetraVector& operator= (const EpetraVector& x)
    
        Assignment operator

    .. cpp:function:: explicit EpetraVector(boost::shared_ptr<Epetra_FEVector> vector)
    
        Create vector view from given Epetra_FEVector pointer

    .. cpp:function:: explicit EpetraVector(const Epetra_Map& map)
    
        Create vector from given Epetra_Map

    .. cpp:function:: explicit EpetraVector(uint N, std::string type="global")
    
        Create vector of size N

    .. cpp:function:: virtual EpetraVector* copy() const
    
        Return copy of tensor

    .. cpp:function:: virtual LinearAlgebraFactory& factory() const
    
        Return linear algebra backend factory

    .. cpp:function:: virtual const EpetraVector& operator*= (const GenericVector& x)
    
        Multiply vector by another vector pointwise

    .. cpp:function:: virtual const EpetraVector& operator*= (double a)
    
        Multiply vector by given number

    .. cpp:function:: virtual const EpetraVector& operator+= (const GenericVector& x)
    
        Add given vector

    .. cpp:function:: virtual const EpetraVector& operator-= (const GenericVector& x)
    
        Subtract given vector

    .. cpp:function:: virtual const EpetraVector& operator/= (double a)
    
        Divide vector by given number

    .. cpp:function:: virtual const EpetraVector& operator= (const GenericVector& x)
    
        Assignment operator

    .. cpp:function:: virtual const EpetraVector& operator= (double a)
    
        Assignment operator

    .. cpp:function:: virtual double inner(const GenericVector& vector) const
    
        Return inner product with given vector

    .. cpp:function:: virtual double max() const
    
        Return maximum value of vector

    .. cpp:function:: virtual double min() const
    
        Return minimum value of vector

    .. cpp:function:: virtual double norm(std::string norm_type) const
    
        Return norm of vector

    .. cpp:function:: virtual double sum() const
    
        Return sum of values of vector

    .. cpp:function:: virtual double sum(const Array<uint>& rows) const
    
        Return sum of selected rows in vector

    .. cpp:function:: virtual std::pair<uint, uint> local_range() const
    
        Return local ownership range of a vector

    .. cpp:function:: virtual std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: virtual uint size() const
    
        Return size of vector

    .. cpp:function:: virtual void add(const double* block, uint m, const uint* rows)
    
        Add block of values

    .. cpp:function:: virtual void add_local(const Array<double>& values)
    
        Add all values to each entry on local process

    .. cpp:function:: virtual void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: virtual void axpy(double a, const GenericVector& x)
    
        Add multiple of given vector (AXPY operation)

    .. cpp:function:: virtual void gather(GenericVector& x, const Array<uint>& indices) const
    
        Gather entries into local vector x

    .. cpp:function:: virtual void get(double* block, uint m, const uint* rows) const
    
        Get block of values

    .. cpp:function:: virtual void get_local(Array<double>& values) const
    
        Get all values on local process

    .. cpp:function:: virtual void resize(uint N)
    
        Resize vector to size N

    .. cpp:function:: virtual void set(const double* block, uint m, const uint* rows)
    
        Set block of values

    .. cpp:function:: virtual void set_local(const Array<double>& values)
    
        Set all values on local process

    .. cpp:function:: virtual void zero()
    
        Set all entries to zero and keep any sparse structure

    .. cpp:function:: virtual ~EpetraVector()
    
        Destructor

    .. cpp:function:: void reset(const Epetra_Map& map)
    
        Reset Epetra_FEVector

