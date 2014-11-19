
.. Documentation for the header file dolfin/log/Table.h

.. _programmers_reference_cpp_log_table:

Table.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Table

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class provides storage and pretty-printing for tables.
    Example usage:
    
      Table table("Timings");
    
      table("uBLAS",  "Assemble") = 0.010;
      table("uBLAS",  "Solve")    = 0.020;
      table("PETSc",  "Assemble") = 0.011;
      table("PETSc",  "Solve")    = 0.019;
      table("Epetra", "Assemble") = 0.012;
      table("Epetra", "Solve")    = 0.018;
    
      info(table);


    .. cpp:function:: Table(std::string title="", bool right_justify=true)
    
        Create empty table


    .. cpp:function:: TableEntry operator() (std::string row, std::string col)
    
        Return table entry


    .. cpp:function:: void set(std::string row, std::string col, int value)
    
        Set value of table entry


    .. cpp:function:: void set(std::string row, std::string col, std::size_t value)
    
        Set value of table entry


    .. cpp:function:: void set(std::string row, std::string col, double value)
    
        Set value of table entry


    .. cpp:function:: void set(std::string row, std::string col, std::string value)
    
        Set value of table entry


    .. cpp:function:: std::string get(std::string row, std::string col) const
    
        Get value of table entry


    .. cpp:function:: double get_value(std::string row, std::string col) const
    
        Get value of table entry


    .. cpp:function:: std::string title() const
    
        Return table title


    .. cpp:function:: const Table& operator= (const Table& table)
    
        Assignment operator


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: std::string str_latex() const
    
        Return informal string representation for LaTeX


.. cpp:class:: TableEntry

    This class represents an entry in a Table


    .. cpp:function:: TableEntry(std::string row, std::string col, Table& table)
    
        Create table entry


    .. cpp:function:: const TableEntry& operator= (std::size_t value)
    
        Assign value to table entry


    .. cpp:function:: const TableEntry& operator= (int value)
    
        Assign value to table entry


    .. cpp:function:: const TableEntry& operator= (double value)
    
        Assign value to table entry


    .. cpp:function:: const TableEntry& operator= (std::string value)
    
        Assign value to table entry


    .. cpp:function:: operator std::string() const
    
        Cast to entry value


