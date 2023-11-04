#pragma once
#include "SupplyCeil.h"
#include <vector>
#include <string>

class TransportTask
{
public:
    // Constructor for the TransportTask class, which initializes its members.
    TransportTask() = default;

    // Delete the copy constructor and copy assignment operator to prevent copying.
    TransportTask(const TransportTask&) = delete;
    TransportTask& operator=(const TransportTask&) = delete;

    // Public methods:

    // Method to generate a transportation problem instance from an input file.
    bool generate(const std::string& input);

    // Method to run the transportation problem solver.
    void run();

    // Method to print the result of the transportation problem.
    void printResult(const std::string& output = "") const;

    // Constant static member representing an empty cell in the transportation matrix.
    inline static const SupplyCeil NOT = {};

private:
    // Private methods:

    // Method to apply the Minimum Cost Method for finding an initial feasible solution.
    void minimumCostMethod();

    // Method to convert the transportation matrix to a vector of SupplyCeil elements.
    std::vector<SupplyCeil> matrixToVector() const;

    // Method to find neighboring cells to a given cell in the matrix.
    std::vector<SupplyCeil> neighbours(const SupplyCeil& sc, const std::vector<SupplyCeil>& vsc) const;

    // Method to compute a cycle path for a given cell.
    std::vector<SupplyCeil> getCiclePath(const SupplyCeil& s) const;

    // Method to check if the current transportation plan is correct.
    bool checkCorrectPlan() const;

    // Method to fix a degenerate transportation plan if needed.
    void fixDegeneratePlan();

    // Method to find the entering cell using the potential method.
    SupplyCeil potentialMethod() const;

    // Method to resolve the transportation plan iteratively.
    void resolvePlan();

    // Private data members:

    std::vector<int> offer, demand;  // Supplies (A) and demands (B)
    std::vector<std::vector<double>> costs;  // Transportation costs
    std::vector<std::vector<SupplyCeil>> transportMatrix;  // Transportation matrix
    size_t iterations = 0;  // Iteration count
};
