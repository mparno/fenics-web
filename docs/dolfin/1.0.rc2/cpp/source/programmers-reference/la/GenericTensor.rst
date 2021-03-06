
.. Documentation for the header file dolfin/la/GenericTensor.h

.. _programmers_reference_cpp_la_generictensor:

GenericTensor.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericTensor

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class defines a common interface for arbitrary rank tensors.


    .. cpp:function:: bool distributed() const = 0
    
        Return true if tensor is distributed


    .. cpp:function:: void init(const GenericSparsityPattern& sparsity_pattern) = 0
    
        Initialize zero tensor using sparsity pattern


    .. cpp:function:: GenericTensor* copy() const = 0
    
        Return copy of tensor


    .. cpp:function:: uint rank() const = 0
    
        Return tensor rank (number of dimensions)


    .. cpp:function:: uint size(uint dim) const = 0
    
        Return size of given dimension


    .. cpp:function:: std::pair<uint, uint> local_range(uint dim) const = 0
    
        Return local ownership range


    .. cpp:function:: void get(double* block, const uint* num_rows, const uint * const * rows) const = 0
    
        Get block of values


    .. cpp:function:: void set(const double* block, const uint* num_rows, const uint * const * rows) = 0
    
        Set block of values


    .. cpp:function:: void add(const double* block, const std::vector<const std::vector<uint>* >& rows) = 0
    
        Add block of values


    .. cpp:function:: void add(const double* block, const std::vector<std::vector<uint> >& rows) = 0
    
        Add block of values


    .. cpp:function:: void add(const double* block, const uint* num_rows, const uint * const * rows) = 0
    
        Add block of values


    .. cpp:function:: void zero() = 0
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode) = 0
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


    .. cpp:function:: LinearAlgebraFactory& factory() const = 0
    
        Return linear algebra backend factory


    .. cpp:function:: const T& down_cast() const
    
        Cast a GenericTensor to its derived class (const version)


    .. cpp:function:: T& down_cast()
    
        Cast a GenericTensor to its derived class (non-const version)


    .. cpp:function:: static boost::shared_ptr<X> down_cast(const boost::shared_ptr<Y> A)
    
        Cast a GenericTensor shared ptr to its derived class. Caller
        must check for success (returns null if cast fails).


    .. cpp:function:: bool has_type() const
    
        Check whether the GenericTensor instance matches a specific type


    .. cpp:function:: const GenericTensor* instance() const
    
        Return concrete instance / unwrap (const version)


    .. cpp:function:: GenericTensor* instance()
    
        Return concrete instance / unwrap (non-const version)


    .. cpp:function:: boost::shared_ptr<const GenericTensor> shared_instance() const
    
        Return concrete shared ptr instance / unwrap (const version)


    .. cpp:function:: boost::shared_ptr<GenericTensor> shared_instance()
    
        Return concrete shared ptr instance / unwrap


    .. cpp:function:: const GenericTensor& operator= (const GenericTensor& x)
    
        Assignment (must be overloaded by subclass)


