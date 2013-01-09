
.. Documentation for the header file dolfin/la/PETScBaseMatrix.h

.. _programmers_reference_cpp_la_petscbasematrix:

PETScBaseMatrix.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScMatrixDeleter

.. cpp:class:: PETScBaseMatrix

    *Parent class(es)*
    
        * :cpp:class:`PETScObject`
        
        * :cpp:class:`Variable`
        
    This class is a base class for matrices that can be used in
    PETScKrylovSolver.


    .. cpp:function:: PETScBaseMatrix()
    
        Constructor


    .. cpp:function:: PETScBaseMatrix(boost::shared_ptr<Mat> A)
    
        Constructor


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return number of rows (dim = 0) or columns (dim = 1)


    .. cpp:function:: std::pair<std::size_t, std::size_t> local_range(std::size_t dim) const
    
        Return local range along dimension dim


    .. cpp:function:: void resize(GenericVector& z, std::size_t dim) const
    
        Resize matrix to be compatible with the matrix-vector product
        y = Ax. In the parallel case, both size and layout are
        important.
        
        *Arguments*
            dim (std::size_t)
                The dimension (axis): dim = 0 --> z = y, dim = 1 --> z = x


    .. cpp:function:: boost::shared_ptr<Mat> mat() const
    
        Return PETSc Mat pointer


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


