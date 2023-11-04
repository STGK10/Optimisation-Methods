#include "TransportTask.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iomanip>
#include <cmath>
#include <sstream>


bool TransportTask::generate(const std::string& input)
{
	// Open the input file stream using the given file path
	std::ifstream file(input);

	// Check if the file is open and can be read
	if (!file.is_open()) {
		return false;
	}

	// Declare variables to store the number of sources and destinations
	size_t numSources, numDest;

	// Read the number of sources and destinations from the file
	file >> numSources >> numDest;

	// Check if either the number of sources or destinations is zero
	if (!numSources || !numDest) {
		return false;
	}

	// Resize the 'offer' and 'demand' vectors to store the supply and demand values
	offer.resize(numSources);
	demand.resize(numDest);

	// Read supply values for each source and demand values for each destination
	for (auto& num : offer) {
		file >> num;
	}
	for (auto& num : demand) {
		file >> num;
	}

	// Calculate the total supply and demand
	size_t cSrc = std::accumulate(offer.cbegin(), offer.cend(), 0);
	size_t cDst = std::accumulate(demand.cbegin(), demand.cend(), 0);

	// If there is an imbalance (total supply doesn't match total demand),
	// introduce a fictitious supplier or consumer to balance the problem
	if (cSrc > cDst) {
		demand.push_back(cSrc - cDst);
	}
	else if (cSrc < cDst) {
		offer.push_back(cDst - cSrc);
	}

	// Loop through the sources and destinations to initialize the costs and transportMatrix
	for (size_t i = 0; i < offer.size(); i++) {
		// Initialize vectors to store cost rates and supply ceils for each destination
		std::vector<double> rates(demand.size(), 0);
		std::vector<SupplyCeil> sc(demand.size());

		// Read cost rates if there are any sources
		if (i < numSources) {
			for (size_t j = 0; j < numDest; j++) {
				file >> rates[j];
			}
		}

		// Store the cost rates and supply ceils for this source
		costs.push_back(rates);
		transportMatrix.push_back(sc);
	}

	// Close the input file
	file.close();

	// Return true to indicate that data initialization was successful
	return true;
}

void TransportTask::run()
{
	//northWestWay();
	minimumCostMethod();
	resolvePlan();
	printResult();
}


void TransportTask::minimumCostMethod() {
	std::vector<std::vector<double>> reducedCosts = costs; // Create a copy of the cost matrix

	while (true) {
		// Find the cell with the minimum cost in the reduced cost matrix
		double minCost = std::numeric_limits<double>::max();
		size_t minRow = 0, minCol = 0;

		for (size_t i = 0; i < offer.size(); i++) {
			for (size_t j = 0; j < demand.size(); j++) {
				if (reducedCosts[i][j] < minCost) {
					minCost = reducedCosts[i][j];
					minRow = i;
					minCol = j;
				}
			}
		}

		// Determine the allocation amount
		int allocation = std::min(offer[minRow], demand[minCol]);

		// Update the transportation matrix with the allocation
		if (allocation) {
			transportMatrix[minRow][minCol] = SupplyCeil(allocation, costs[minRow][minCol], minRow, minCol);
		}
		

		// Update supply and demand
		offer[minRow] -= allocation;
		demand[minCol] -= allocation;

		// Update reduced costs
		reducedCosts[minRow][minCol] = std::numeric_limits<double>::max();

		// Check for completion
		bool allAllocated = true;
		for (size_t i = 0; i < offer.size(); i++) {
			for (size_t j = 0; j < demand.size(); j++) {
				if (reducedCosts[i][j] < std::numeric_limits<double>::max()) {
					allAllocated = false;
					break;
				}
			}
			if (!allAllocated) {
				break;
			}
		}

		if (allAllocated) {
			break;
		}
	}
}



/*
// Выбираем первый опорный план, как и в Симплекс-методе
void TransportTask::northWestWay()
{
	size_t northWest = 0;
	
	for (size_t i = 0; i < offer.size(); i++)
	{
		for (size_t j = northWest; j < demand.size(); j++)
		{
			int resource = std::min(offer[i], demand[j]);
			
			if (resource >= 0)
			{
				transportMatrix[i][j] = SupplyCeil(resource, costs[i][j], i, j);
				offer[i]  -= resource;
				demand[j] -= resource;
				
				// Закрываем столбец
				if (offer[i] == 0)
				{
					northWest = j;
					break;
				}
			}
		}
	}
}

*/
std::vector<SupplyCeil> TransportTask::matrixToVector() const
{
	// Create an empty vector to store the converted supply ceils
	std::vector<SupplyCeil> result;

	// Iterate over each row in the 'transportMatrix'
	for (const auto& row : transportMatrix)
	{
		// Iterate over each element (SupplyCeil) in the current row
		for (const auto& sc : row)
		{
			// Check if the current SupplyCeil is not equal to the special value TransportTask::NOT
			if (sc != TransportTask::NOT) {
				// If it's not marked as NOT, add it to the result vector
				result.push_back(sc);
			}
		}
	}

	// Return the result vector containing the non-NOT SupplyCeil objects
	return result;
}



std::vector<SupplyCeil> TransportTask::neighbours(const SupplyCeil& sc, const std::vector<SupplyCeil>& vsc) const
{
	// Create an empty vector to store the neighboring SupplyCeil objects
	std::vector<SupplyCeil> neighbour(2);

	// Iterate over each SupplyCeil object in the input vector 'vsc'
	for (const auto& supp : vsc)
	{
		// Check if the current SupplyCeil 'supp' is not equal to the reference SupplyCeil 'sc'
		if (supp != sc)
		{
			// Check if the row of 'supp' matches the row of 'sc'
			if (supp.row == sc.row && neighbour[0] == TransportTask::NOT)
			{
				// If the condition is met, set the first neighbor in 'neighbour' to 'supp'
				neighbour[0] = supp;
			}
			// Check if the column of 'supp' matches the column of 'sc'
			else if (supp.col == sc.col && neighbour[1] == TransportTask::NOT)
			{
				// If the condition is met, set the second neighbor in 'neighbour' to 'supp'
				neighbour[1] = supp;
			}

			// Check if both neighbors have been found (i.e., not equal to TransportTask::NOT)
			if (neighbour[0] != TransportTask::NOT && neighbour[1] != TransportTask::NOT)
			{
				// If both neighbors are found, exit the loop
				break;
			}
		}
	}

	// Return the 'neighbour' vector containing the two nearest neighbors (if found)
	return neighbour;
}

//finding and returning a cycle path in a transportation problem based on a reference
std::vector<SupplyCeil> TransportTask::getCiclePath(const SupplyCeil& sc) const
{
	// Create a vector 'path' and initialize it by converting the 'transportMatrix' into a 1D vector
	auto path = matrixToVector();

	// Insert the reference 'sc' at the beginning of the 'path' vector
	path.insert(path.begin(), sc);

	// Start a loop to identify the cycle path
	size_t previous = 0;
	do
	{
		previous = path.size();

		// Use the remove_if algorithm to filter elements in 'path' that do not have neighbors in horizontal or vertical directions
		auto iter = std::remove_if(path.begin(), path.end(), [&path, this](const SupplyCeil& c) {
			auto neighbour = neighbours(c, path);
			return neighbour[0] == TransportTask::NOT || neighbour[1] == TransportTask::NOT;
			});

		// Erase elements that were marked for removal by the remove_if algorithm
		path.erase(iter, path.end());

	} while (previous != path.size());

	// Create a result vector with the same size as 'path' and initialize it with TransportTask::NOT
	std::vector<SupplyCeil> result(path.size(), TransportTask::NOT);

	// Initialize 'prev' with the reference 'sc'
	auto prev = sc;

	// Loop through the result vector to construct the correct cycle path
	for (size_t i = 0; i < result.size(); i++)
	{
		// Set the current position in the result vector to 'prev'
		result[i] = prev;

		// Update 'prev' with the next neighbor from 'path', alternating between neighbors
		prev = neighbours(prev, path)[i % 2];
	}

	// Return the result vector containing the cycle path
	return result;
}

// addressing a degenerate transportation problem by ensuring that the solution plan is non-degenerate
void TransportTask::fixDegeneratePlan()
{
	// Loop through the offer and demand values to inspect the transportation matrix
	for (size_t i = 0; i < offer.size(); i++)
	{
		for (size_t j = 0; j < demand.size(); j++)
		{
			// Check if the current cell in the transportation matrix is marked as TransportTask::NOT
			if (transportMatrix[i][j] == TransportTask::NOT)
			{
				// To fix a degenerate plan, we need to check for the acyclicity property

				// Create a fictional supply ceil with no allocation (0 resources) and associated cost
				SupplyCeil fictional(0, costs[i][j], i, j);

				// Check if the cycle path starting from the fictional supply ceil is empty
				if (getCiclePath(fictional).empty())
				{
					// If there are no cycles in the cycle path, set the transportation matrix cell to the fictional supply ceil
					transportMatrix[i][j] = fictional;

					// Return from the function, as we have resolved the degeneracy
					return;
				}
			}
		}
	}
}

// checking whether the current plan is correct in the context of the transportation problem.
bool TransportTask::checkCorrectPlan() const{
	return offer.size() + demand.size() - 1 == matrixToVector().size();
}

//calculates the potentials (also known as dual variables) for the transportation problem and uses them to find a cell with the maximum positive deviation from the current solution
SupplyCeil TransportTask::potentialMethod() const
{
	// Create vectors 'u' and 'v' to store potentials and initialize them with NaN values
	std::vector<double> u(offer.size(), std::numeric_limits<double>::quiet_NaN());
	std::vector<double> v(demand.size(), std::numeric_limits<double>::quiet_NaN());

	// Initialize the first element of 'u' with 0.0
	u[0] = 0.0;

	// Enter a loop that calculates the potentials until certain conditions are met
	while (true)
	{
		// Check if any element in 'u' or 'v' is NaN (Not-a-Number)
		auto isNanU = std::any_of(u.cbegin(), u.cend(), std::isnan<double>);
		auto isNanV = std::any_of(v.cbegin(), v.cend(), std::isnan<double>);

		// If both 'u' and 'v' have no NaN elements, exit the loop
		if (!isNanU && !isNanV)
		{
			break;
		}

		// Iterate through the offer and demand values to update potentials 'u' and 'v'
		for (size_t i = 0; i < offer.size(); i++)
		{
			for (size_t j = 0; j < demand.size(); j++)
			{
				// Check if the current cell in the transportation matrix is marked as TransportTask::NOT
				if (transportMatrix[i][j] == TransportTask::NOT) {
					continue;
				}

				// Check if either 'u[i]' or 'v[j]' is not NaN
				if (std::isnan(u[i]) && std::isnan(v[j])) {
					continue;
				}

				// Update 'u[i]' if it's NaN based on the value of 'v[j]' and the cost matrix
				if (std::isnan(u[i])) {
					u[i] = v[j] - costs[i][j];
				}
				// Update 'v[j]' if it's NaN based on the value of 'u[i]' and the cost matrix
				else if (std::isnan(v[j])) {
					v[j] = costs[i][j] + u[i];
				}
			}
		}
	}

	// Initialize variables to find the cell with maximum positive deviation
	SupplyCeil result = TransportTask::NOT;
	double a_ij = std::numeric_limits<double>::min();

	// Iterate through the offer and demand values to find the cell with maximum positive deviation
	for (size_t i = 0; i < offer.size(); i++)
	{
		for (size_t j = 0; j < demand.size(); j++)
		{
			// Check if the current cell in the transportation matrix is not marked as TransportTask::NOT
			if (transportMatrix[i][j] != TransportTask::NOT) {
				continue;
			}

			// Calculate the difference between 'v[j]', 'u[i]', and the cost matrix value 'costs[i][j]'
			double difference = v[j] - u[i] - costs[i][j];

			// Choose the cell that corresponds to the maximum positive difference
			if (difference > 0.0)
			{
				if (difference > a_ij)
				{
					a_ij = difference;
					result = SupplyCeil(0, costs[i][j], i, j);
				}
			}
		}
	}

	// Return the selected cell with the maximum positive deviation
	return result;
}

/*
* Комментарии по поводу схожести странспортного алгоритма и симплекс-метода.
* 
* В симплекс методе все переменные делились на два непересекающихся подмножетсва: базисных и свободных переменных
* Для оптимизации опорного плана при определенных условиях мы меняли местами переменные из различных множеств.
* Одна базисная переменная становилась свободной, а свободная - базисной, после чего значение целевой функции пересчитывалось, алгоритм переходил к следующему шагу.
* 
* Аналогичным образом мы поступаем в транспортной задаче. Чтобы подобрать оптимальный план, мы последовательно меняем базисные и свободные переменные.
* Подбор таких переменных осуществляется посредством метода потенциалов, а также анализа "циклов".
*/

//implements the process of finding an optimal solution to a transportation problem using the modified distribution method
void TransportTask::resolvePlan()
{
	// The while loop continues until the termination condition is met.
	while (true)
	{
		// Check if the current plan is correct (whether it satisfies certain conditions)
		if (!checkCorrectPlan()) {
			// If the plan is degenerate, fix it to make it non-degenerate
			fixDegeneratePlan();
		}

		// Find a cell with the maximum positive deviation (potential method)
		SupplyCeil leavingCeil = potentialMethod();

		// If no cell with a positive deviation is found, the algorithm terminates
		if (leavingCeil == TransportTask::NOT) {
			return;
		}

		// Find a cycle path in the transportation matrix starting from the leaving cell
		auto path = getCiclePath(leavingCeil);

		// Initialize variables to find the minimum resource value in the cycle path
		double minRes = std::numeric_limits<double>::max();
		bool isPlus = true;

		// Find the minimum resource value (minRes) in the cycle path based on the "plus" and "minus" rule
		for (const auto& sc : path)
		{
			if (!isPlus)
			{
				// Check if the current cell is not the leaving cell and has a smaller resource value
				if (sc != leavingCeil && transportMatrix[sc.row][sc.col].resources < minRes) {
					minRes = sc.resources;
				}
			}

			isPlus = !isPlus;
		}

		// Update the leaving cell with the new resource value
		transportMatrix[leavingCeil.row][leavingCeil.col] = SupplyCeil(minRes, costs[leavingCeil.row][leavingCeil.col], leavingCeil.row, leavingCeil.col);

		// Perform a cyclic recalculation by updating the resources of cells in the cycle path
		// and removing cells with nearly zero resources
		bool isTrash = true;
		isPlus = false;
		for (const auto& sc : path)
		{
			isPlus = !isPlus;
			if (sc.row == leavingCeil.row && sc.col == leavingCeil.col) {
				continue;
			}

			// Update the resource value based on the "plus" and "minus" rule
			transportMatrix[sc.row][sc.col].resources += (2.0 * isPlus - 1.0) * minRes;

			// If the updated resource value is very close to zero, mark the cell as TransportTask::NOT
			if (transportMatrix[sc.row][sc.col].resources < std::numeric_limits<double>::epsilon() && isTrash) {
				transportMatrix[sc.row][sc.col] = TransportTask::NOT;
				isTrash = false;
			}
		}

		// Increment the number of iterations
		iterations++;
	}
}


// printing the result of the transportation problem, including the optimized plan and total cost.
void TransportTask::printResult(const std::string& output) const
{
	// Create an output string stream to build the result output
	std::ostringstream oss;

	// Add a header to the result output, indicating the number of iterations
	oss << "Optimal plan was goal with " << iterations << ": \n\n";
	double totalCost = 0.0;  // Initialize the total cost variable to zero

	// Loop through the rows (sources)
	for (size_t i = 0; i < offer.size(); i++)
	{
		oss << std::setw(4) << "A" << (i + 1);  // Label each row as "A1," "A2," etc.

		// Loop through the columns (destinations)
		for (size_t j = 0; j < demand.size(); j++)
		{
			const auto& sc = transportMatrix[i][j];  // Get the supply ceil for the current cell

			// Check if the supply ceil is not marked as TransportTask::NOT and belongs to the current row and column
			if (sc != TransportTask::NOT && sc.row == i && sc.col == j)
			{
				oss << std::setw(4) << sc.resources;  // Print the resource value in the cell
				totalCost += sc.resources * sc.costPerTrans;  // Update the total cost
			}
			else
			{
				oss << std::setw(4) << "-";  // Print a dash "-" if the cell is not part of the solution
			}
		}

		oss << "\n";  // Move to the next line after processing all columns for the current row
	}

	// Print column labels ("B1," "B2," etc.) in the result output
	oss << std::setw(3) << "     ";
	for (size_t j = 0; j < demand.size(); j++) {
		oss << std::setw(3) << "B" << (j + 1);  // Label each column as "B1," "B2," etc.
	}

	oss << "\n\nTotal cost: " << totalCost << "\n";  // Print the total cost at the end of the output

	if (!output.empty())
	{
		// If an output file name is provided, write the result to that file
		std::ofstream ofs(output);
		ofs << oss.str();
		ofs.close();
	}
	else
	{
		// If no output file name is provided, print the result to the console
		std::cout << oss.str() << std::flush;
	}
}
