
.. Documentation for the header file dolfin/la/VectorSpaceBasis.h

.. _programmers_reference_cpp_la_vectorspacebasis:

VectorSpaceBasis.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: VectorSpaceBasis

    This class defines a basis for vector spaces,
    typically used for expressing nullspaces, transpose nullspaces
    and near nullspaces of singular operators


    .. cpp:function:: VectorSpaceBasis(const std::vector<std::shared_ptr<GenericVector> > basis)
    
        Constructor


    .. cpp:function:: bool is_orthonormal() const
    
        Test if basis is orthonormal


    .. cpp:function:: bool is_orthogonal() const
    
        Test if basis is orthogonal


    .. cpp:function:: void orthogonalize(GenericVector& x) const
    
        Orthogonalize x with respect to basis


    .. cpp:function:: std::size_t dim() const
    
        Dimension of the basis


    .. cpp:function:: std::shared_ptr<const GenericVector> operator[] (std::size_t i) const
    
        Get a particular basis vector


