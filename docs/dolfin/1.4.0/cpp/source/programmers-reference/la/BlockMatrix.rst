
.. Documentation for the header file dolfin/la/BlockMatrix.h

.. _programmers_reference_cpp_la_blockmatrix:

BlockMatrix.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BlockMatrix

    .. cpp:function:: void set_block(std::size_t i, std::size_t j, std::shared_ptr<GenericMatrix> m)
    
        Set block


    .. cpp:function:: std::shared_ptr<const GenericMatrix> get_block(std::size_t i, std::size_t j) const
    
        Get block (const version)


    .. cpp:function:: std::shared_ptr<GenericMatrix> get_block(std::size_t i, std::size_t j)
    
        Get block


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: void zero()
    
        Set all entries to zero and keep any sparse structure


    .. cpp:function:: void apply(std::string mode)
    
        Finalize assembly of tensor


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void mult(const BlockVector& x, BlockVector& y, bool transposed=false) const
    
        Matrix-vector product, y = Ax


    .. cpp:function:: std::shared_ptr<GenericMatrix> schur_approximation(bool symmetry=true) const
    
        Create a crude explicit Schur approximation  of S = D - C A^-1
        B of  (A B; C  D) If symmetry !=  0, then the  caller promises
        that B = symmetry * transpose(C).


