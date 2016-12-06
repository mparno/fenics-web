
.. Documentation for the header file dolfin/la/test_nullspace.h

.. _programmers_reference_cpp_la_test_nullspace:

test_nullspace.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: bool in_nullspace(const GenericLinearOperator& A, const VectorSpaceBasis& x, std::string type="right")

    Check whether a vector space basis is in the nullspace of a
    given operator. The string option 'type' can be "right" for the
    right nullspace (Ax=0) or "left" for the left nullspace (A^Tx =
    0). To test the left nullspace, A must also be of type
    GenericMatrix.


