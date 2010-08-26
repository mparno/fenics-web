.. Documentation for the header file dolfin/mesh/MeshFunction.h

.. _programmers_reference_cpp_mesh_Mesh:

MeshFunction.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`Variable`
        
        A MeshFunction is a function that can be evaluated at a set of
        mesh entities. A MeshFunction is discrete and is only defined
        at the set of mesh entities of a fixed topological dimension.
        A MeshFunction may for example be used to store a global
        numbering scheme for the entities of a (parallel) mesh, marking
        sub domains or boolean markers for mesh refinement.

    .. cpp:function:: MeshFunction() :
                       Variable("f", "unnamed MeshFunction"),
                       _values(0), _mesh(0), _dim(0), _size(0)
    
        Create empty mesh function

    .. cpp:function:: MeshFunction(const Mesh& mesh) :
                                   Variable("f", "unnamed MeshFunction"),
                                   _values(0), _mesh(&mesh), _dim(0), _size(0)
    
        Create empty mesh function on given mesh

    .. cpp:function:: MeshFunction(const Mesh& mesh, const std::string filename) :
                                   Variable("f", "unnamed MeshFunction"),
                                   _values(0), _mesh(&mesh), _dim(0), _size(0)
    
        Create function from data file

    .. cpp:function:: MeshFunction(const Mesh& mesh, uint dim) :
                                   Variable("f", "unnamed MeshFunction"),
                                   _values(0), _mesh(&mesh), _dim(0), _size(0)
    
        Create mesh function on given mesh of given dimension

    .. cpp:function:: MeshFunction(const Mesh& mesh, uint dim, const T& value) :
                                   Variable("f", "unnamed MeshFunction"),
                                   _values(0), _mesh(&mesh), _dim(0), _size(0)
    
        Create mesh function on given mesh of given dimension and initialise
        to a value

    .. cpp:function:: MeshFunction(const MeshFunction<T>& f) :
                                   Variable("f", "unnamed MeshFunction"),
                                   _values(0), _mesh(0), _dim(0), _size(0)
    
        Copy constructor

    .. cpp:function:: T& operator[] (const MeshEntity& entity)
    
        Return value at given entity

    .. cpp:function:: T& operator[] (uint index)
    
        Return value at given index

    .. cpp:function:: T* _values
    
        Values at the set of mesh entities

    .. cpp:function:: T* values()
    
        Return array of values

    .. cpp:function:: const Mesh& mesh() const
    
        Return mesh associated with mesh function

    .. cpp:function:: const Mesh* _mesh
    
        The mesh

    .. cpp:function:: const MeshFunction<T>& operator= (const MeshFunction<T>& f)
    
        Assign mesh function

    .. cpp:function:: const MeshFunction<T>& operator= (const T& value)
    
        Set all values to given value

    .. cpp:function:: const T& operator[] (const MeshEntity& entity) const
    
        Return value at given entity (const version)

    .. cpp:function:: const T& operator[] (uint index) const
    
        Return value at given index  (const version)

    .. cpp:function:: const T* values() const
    
        Return array of values

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: uint _dim
    
        Topological dimension

    .. cpp:function:: uint _size
    
        Number of mesh entities

    .. cpp:function:: uint dim() const
    
        Return topological dimension

    .. cpp:function:: uint size() const
    
        Return size (number of entities)

    .. cpp:function:: void init(const Mesh& mesh, uint dim)
    
        Initialize mesh function for given topological dimension

    .. cpp:function:: void init(const Mesh& mesh, uint dim, uint size)
    
        Initialize mesh function for given topological dimension of given size

    .. cpp:function:: void init(uint dim)
    
        Initialize mesh function for given topological dimension

    .. cpp:function:: void init(uint dim, uint size)
    
        Initialize mesh function for given topological dimension of given size

    .. cpp:function:: void set_all(const T& value)
    
        Set all values to given value

    .. cpp:function:: ~MeshFunction()
    
        Destructor

