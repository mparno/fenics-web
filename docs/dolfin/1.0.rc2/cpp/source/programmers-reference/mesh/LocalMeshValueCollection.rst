
.. Documentation for the header file dolfin/mesh/LocalMeshValueCollection.h

.. _programmers_reference_cpp_mesh_localmeshvaluecollection:

LocalMeshValueCollection.h
==========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalMeshValueCollection

    This class stores mesh data on a local processor corresponding
    to a portion of a MeshValueCollection.


    .. cpp:function:: LocalMeshValueCollection(const MeshValueCollection<T>& values, uint dim)
    
        Create local mesh data for given LocalMeshValueCollection


    .. cpp:function:: uint dim () const
    
        Return dimension of cell entity


    .. cpp:function:: const std::vector<std::pair<std::pair<uint, uint>, T> >& values() const
    
        Return data


