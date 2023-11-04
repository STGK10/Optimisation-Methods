#pragma once
#include <cstddef>


//represents a cell or ceil in a transportation problem, storing information about the amount of resources, the cost per unit of transportation, and the row and column indices of the cell
struct SupplyCeil
{
    double resources;        // The amount of resources in the ceil.
    double costPerTrans;    // The cost per unit of transportation in the ceil.
    size_t row, col;        // The row and column indices of the ceil.

    // Default constructor for SupplyCeil
    SupplyCeil() : resources(0), costPerTrans(0), row(~0), col(~0) {}

    // Parameterized constructor for SupplyCeil
    SupplyCeil(double res, double cpt, size_t r, size_t c)
        : resources(res), costPerTrans(cpt), row(r), col(c) {}

    // Custom equality operator for SupplyCeil
    friend bool operator==(const SupplyCeil& lhs, const SupplyCeil& rhs)
    {
        return lhs.resources == rhs.resources &&         // Check if resources are equal
            lhs.costPerTrans == rhs.costPerTrans && // Check if costPerTrans is equal
            lhs.row == rhs.row &&                     // Check if row indices are equal
            lhs.col == rhs.col;                      // Check if column indices are equal
    }

    // Custom inequality operator for SupplyCeil
    friend bool operator!=(const SupplyCeil& lhs, const SupplyCeil& rhs)
	{
    	return !(lhs == rhs);
	}
};
