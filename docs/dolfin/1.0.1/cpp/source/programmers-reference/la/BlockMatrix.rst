
.. Documentation for the header file dolfin/la/BlockMatrix.h

.. _programmers_reference_cpp_la_blockmatrix:

BlockMatrix.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BlockMatrix

    .. cpp:function:: void set_block(uint i, uint j, boost::shared_ptr<GenericMatrix> m)
    
        Set block


    .. cpp:function:: const boost::shared_ptr<GenericMatrix> get_block(uint i, uint j) const
    
        Get block (const version)


    .. cpp:function:: boost::shared_ptr<GenericMatrix> get_block(uint i, uint j)
    
        Get block


    .. cpp:function:: uint size(uint dim) const
    
        Return size of given dimension


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void mult(const BlockVector& x, BlockVector& y, bool transposed=false) const
    
        Matrix-vector product, y = Ax


    .. cpp:function:: boost::shared_ptr<GenericMatrix> schur_approximation(bool symmetry=true) const
    
        Create a crude explicit Schur approximation of S = D - C A^-1 B of (A B; C D)
        If symmetry != 0, then the caller promises that B = symmetry * transpose(C).


