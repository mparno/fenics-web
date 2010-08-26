.. Documentation for the header file dolfin/la/BlockMatrix.h

.. _programmers_reference_cpp_la_Mesh:

BlockMatrix.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: BlockMatrix

    .. cpp:function:: Matrix& get(uint i, uint j)
    
        Get block

    .. cpp:function:: SubMatrix operator() (uint i, uint j)
    
        Return SubMatrix reference number (i,j)

    .. cpp:function:: const Matrix& get(uint i, uint j) const
    
        Get block (const version)

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: uint size(uint dim) const
    
        Return size of given dimension

    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor

    .. cpp:function:: void mult(const BlockVector& x, BlockVector& y, bool transposed=false) const
    
        Matrix-vector product, y = Ax

    .. cpp:function:: void set(uint i, uint j, Matrix& m)
    
        Set block

    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure

.. cpp:class:: SubMatrix

    .. cpp:function:: const SubMatrix& operator= (Matrix& m)
    
        Assign Matrix to SubMatrix

