#include "LPProblem.h"

void LPProblem::vector_print(vector<double>& v) {
	for (auto& elem : v) {
		cout << elem << " ";
	}
	cout << "\n";
}

void LPProblem::matrix_print(vector<vector<double>>& a) {
	for (int i = 0; i < a.size(); i++) {
		for (auto& elem : a[i]) {
			if (elem >= 0) cout << " ";
			cout << elem << " ";
		}
		cout << "\n";
	}
	cout << "\n";
}

void LPProblem::print() {
	for (int i = 0; i < var_num; i++) {
		if (i != 0) cout << " + ";
		cout << fun_coeffs[i] << "x_" << i + 1;
	}
	cout << " --> " << (cond ? "max\n\n" : "min\n\n");

	for (int i = 0; i < res_num; i++) {
		for (int j = 0; j < var_num; j++) {
			if (j != 0) cout << " + ";
			cout << res_matrix[i][j] << "x_" << j + 1;
		}

		if (signs[i] == 0) cout << " = ";
		else if (signs[i] == 1) cout << " <= ";
		else if (signs[i] == 2) cout << " >= ";

		cout << res_vector[i] << endl;
	}

	for (int i = 0; i < var_num; i++) {
		if (res_signs[i] == 1) cout << "x_" << i + 1 << " >= 0; ";
		if (res_signs[i] == 2) cout << "x_" << i + 1 << " <= 0; ";
	}
	cout << "\n\n";
}

void LPProblem::to_canonical_form() {
	if (cond == 1) {
		for (int i = 0; i < var_num; i++)
			fun_coeffs[i] *= -1;
		cond = 0;
		condChange = true;
	}

	for (int i = 0; i < res_num; i++) {
		if (res_vector[i] < 0) {
			res_vector[i] *= -1;
			for (int j = 0; j < var_num; j++) {
				res_matrix[i][j] *= -1;
			}
			if (signs[i] == 1) signs[i] = 2;
			else if (signs[i] == 2) signs[i] = 1;
		}
	}

	for (int i = 0; i < var_num; i++) {
		if (res_signs[i] == 2) {
			res_signs[i] = 1;
			fun_coeffs[i] *= -1;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j][i] *= -1;
			}
		}
		else if (res_signs[i] == 0) {
			res_signs[i] = 1;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j].push_back(res_matrix[j][var_num - 1]);
				for (int k = var_num - 1; k > i + 1; k--) {
					res_matrix[j][k] = res_matrix[j][k - 1];
				}
			}
			res_signs.push_back(res_signs[var_num - 1]);
			for (int j = var_num - 1; j > i + 1; j--) {
				res_signs[j] = res_signs[j - 1];
			}
			res_signs[i + 1] = 1;
			res_signs[i] = 1;

			fun_coeffs.push_back(fun_coeffs[var_num - 1]);
			for (int j = var_num - 1; j > i + 1; j--) {
				fun_coeffs[j] = fun_coeffs[j - 1];
			}
			fun_coeffs[i + 1] = (-1) * fun_coeffs[i];

			var_num++;
			for (int j = 0; j < res_num; j++) {
				res_matrix[j][i + 1] = (-1) * res_matrix[j][i];
			}
		}

	}


	for (int i = 0; i < res_num; i++) {
		if (signs[i]) {
			fun_coeffs.push_back(0);
			for (int j = 0; j < res_num; j++) {
				if (i == j) res_matrix[j].push_back((signs[i] == 1 ? 1 : -1));
				else res_matrix[j].push_back(0);
			}
			signs[i] = 0;
			res_signs.push_back(1);
			var_num++;
		}
	}
}

vector<double> LPProblem::simplex_method2(vector<double>& sol) {
	if (sol.size() != var_num) exit(-2);
	const double zero = 0.000001;

	// Separate variables into non-basic (Nk) and basic (Lk) variables.
	vector<int> Nk, Lk;
	for (int i = 0; i < var_num; i++) {
		if (sol[i] > zero) Nk.push_back(i);
		else Lk.push_back(i);
	}

	while (true) {
		bool isUnpl = true;
		vector<int> zeroCfs;

		if (Nk.size() > res_num) {
			cout << "Unbounded problem.\n";
			exit(-3);
		}

		// Ensure that there are enough basic variables.
		while (Nk.size() < res_num) {
			isUnpl = false;
			zeroCfs.push_back(Nk.size());
			Nk.push_back(Lk[Lk.size() - 1]);
			Lk.pop_back();
		}

		// Prepare matrices ANk and ALk.
		vector<vector<double>> ANk(res_num);
		vector<vector<double>> ALk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int& cf : Nk) {
				ANk[i].push_back(res_matrix[i][cf]);
			}
			for (int& cf : Lk) {
				ALk[i].push_back(res_matrix[i][cf]);
			}
		}

		// Compute the inverse matrix B.
		vector<vector<double>> B = inverse_matrix(ANk);
		if (B.size() == 1) {
			// If B size is 1, the problem is unbounded.
			Lk.push_back(Nk[Nk.size() - 1]);
			Nk.pop_back();
			Nk.push_back(Lk[0]);
			Lk.erase(Lk.begin());
			continue;
		}

		// Prepare coefficient vectors for non-basic (cNk) and basic (cLk) variables.
		vector<double> cLk;
		vector<double> cNk;
		for (int& cf : Lk) {
			cLk.push_back(fun_coeffs[cf]);
		}
		for (int& cf : Nk) {
			cNk.push_back(fun_coeffs[cf]);
		}

		// Calculate cNkBNk and cNkBNkALk.
		vector<double> cNkBNk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				cNkBNk[i] += cNk[j] * B[j][i];
			}
		}

		vector<double> cNkBNkALk(var_num - res_num);
		for (int i = 0; i < var_num - res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				cNkBNkALk[i] += cNkBNk[j] * ALk[j][i];
			}
		}

		// Calculate dLk (difference between cLk and cNkBNkALk).
		vector<double> dLk(var_num - res_num);
		for (int i = 0; i < var_num - res_num; i++) {
			dLk[i] = cLk[i] - cNkBNkALk[i];
		}

		// Check if the current solution is optimal.
		bool isOpt = true;
		int jk = 0;
		for (int i = 0; i < var_num - res_num; i++) {
			if (dLk[i] <= -zero) {
				isOpt = false;
				jk = i;
			}
		}

		if (isOpt) {
			// If the solution is optimal, set negative values to zero and return the solution.
			for (int i = 0; i < var_num; i++) {
				if (sol[i] < zero) sol[i] = 0;
			}
			return sol;
		}

		// Calculate vector 'uk' and check for infeasibility.
		vector<double> uk(res_num);
		for (int i = 0; i < res_num; i++) {
			for (int j = 0; j < res_num; j++) {
				uk[i] += B[i][j] * ALk[j][jk];
			}
		}
		bool isInf = true;
		vector<int> ik;
		for (int i = 0; i < res_num; i++) {
			if (uk[i] > 0) {
				isInf = false;
				ik.push_back(i);
			}
		}

		if (isInf) {
			cout << "Not a bounded task";
			exit(-4);
		}

		if (isUnpl) {
			// Apply the unplanned movement.
			double theta = sol[Nk[ik[0]]] / uk[ik[0]];
			for (int i = 1; i < ik.size(); i++) {
				double temp = sol[Nk[ik[i]]] / uk[ik[i]];
				if (temp < theta && temp > zero) theta = temp;
			}

			for (int i = 0; i < res_num; i++) {
				sol[Nk[i]] -= theta * uk[i];
			}
			sol[Lk[jk]] += theta;

			Lk.clear();
			Nk.clear();
			for (int i = 0; i < var_num; i++) {
				if (sol[i] > zero) Nk.push_back(i);
				else Lk.push_back(i);
			}
		}
		else {
			// Apply the planned movement.
			bool isNeg = true;
			for (int& cf : zeroCfs) {
				if (uk[cf] > zero) isNeg = false;
			}
			if (isNeg) {
				double theta = INT_MAX;
				for (int i = 0; i < ik.size(); i++) {
					double temp = sol[Nk[ik[i]]] / uk[ik[i]];
					if (temp < theta && temp > zero) theta = temp;
				}

				for (int i = 0; i < res_num; i++) {
					sol[Nk[i]] -= theta * uk[i];
				}
				sol[Lk[jk]] += theta;
				Lk.clear();
				Nk.clear();
				for (int i = 0; i < var_num; i++) {
					if (sol[i] > zero) Nk.push_back(i);
					else Lk.push_back(i);
				}
			}
			else {
				Lk.push_back(Nk[Nk.size() - 1]);
				Nk.pop_back();
				Nk.push_back(Lk[0]);
				Lk.erase(Lk.begin());
			}
		}
	}
}



vector<double> LPProblem::simplex() {
	// Step 1: Create a modified matrix A with additional columns for slack variables.
	vector<vector<double>> A = res_matrix;
	for (int i = 0; i < res_num; i++) {
		for (int j = 0; j < res_num; j++) {
			if (i == j) A[i].push_back(1);
			else A[i].push_back(0);
		}
	}

	// Step 2: Prepare the objective function coefficients vector 'c'.
	vector<double> c;
	for (int i = 0; i < var_num; i++) {
		c.push_back(0);  // Initialize non-basic variables with zero coefficients.
	}
	for (int i = 0; i < res_num; i++) {
		c.push_back(1);  // Initialize basic variables (slack variables) with coefficients 1.
	}

	// Step 3: Prepare the right-hand side vector 'r_vec'.
	vector<double> r_vec = res_vector;

	// Step 4: Prepare the vectors 's' and 'r_s'.
	vector<int> s = { 1, 1 };  // This appears to be related to some specific requirements.
	vector<int> r_s = { 1, 1 };  // Similar to 's', it's not clear from the code.

	// Step 5: Create a new LPProblem instance with the modified objective function and constraints.
	LPProblem* new_lp = new LPProblem(c, A, r_vec, s, r_s, 0);

	// Step 6: Initialize an accumulated solution vector 'acc_sol' with zeros.
	vector<double> acc_sol;
	for (int i = 0; i < var_num; i++) {
		acc_sol.push_back(0);
	}
	for (int i = 0; i < res_num; i++) {
		acc_sol.push_back(res_vector[i]);  // Add the initial values for slack variables.
	}

	// Step 7: Solve the LPProblem using the simplex method and get an initial solution.
	vector<double> solution = new_lp->simplex_method2(acc_sol);

	// Step 8: Remove the values of slack variables from the solution.
	for (int i = 0; i < res_num; i++) {
		solution.pop_back();
	}

	// Step 9: Solve the LPProblem again using the simplex method with the initial solution.
	// This will likely refine the solution to the original problem.
	vector<double> sol = this->simplex_method2(solution);

	// Step 10: Return the final optimal solution to the linear programming problem.
	return sol;
}

LPProblem* LPProblem::dual() {
	int new_var_num = res_num;
	int new_res_num = var_num;
	vector<double> new_fun_coeffs = res_vector;
	vector<double> new_res_vector = fun_coeffs;
	vector<int> new_signs;
	for (auto& sign : res_signs) {
		if (sign == 0) new_signs.push_back(0);
		else if (sign == 2) new_signs.push_back(2);
		else if (sign == 1) new_signs.push_back(1);
	}
	vector<int> new_res_signs;
	for (auto& sign : signs) {
		if (sign == 0) new_res_signs.push_back(0);
		else if (sign == 1) new_res_signs.push_back(2);
		else if (sign == 2) new_res_signs.push_back(1);
	}
	vector<vector<double>> new_res_matrix(new_res_num);
	for (int i = 0; i < new_res_num; i++) {
		for (int j = 0; j < new_var_num; j++) {
			new_res_matrix[i].push_back(res_matrix[j][i]);
		}
	}
	int new_cond = (cond ? 0 : 1);
	return new LPProblem(new_fun_coeffs, new_res_matrix, new_res_vector, new_signs, new_res_signs, new_cond);
}

vector<double> LPProblem::extreme_point_enum() {
	// Step 1: Create a bitmask representing the selection of variables.
	string bitmask(res_num, 1);  // Initialize with '1's (basic variables).
	bitmask.resize(var_num, 0);  // Set the remaining '0's for non-basic variables.

	// Step 2: Initialize a matrix and vectors to store solutions and best solutions.
	vector<vector<double>> matrix(res_num);  // Matrix for constraint coefficients.
	vector<double> solution(res_num);  // Current solution vector.
	vector<double> best_solution(res_num);  // Best solution found so far.
	vector<int> vars(res_num);  // Variables corresponding to the solution.
	vector<int> vars_best(res_num);  // Variables corresponding to the best solution.
	double min_sum = INT_MAX;  // Initialize minimum sum with a large value.

	// Step 3: Enumerate through all possible combinations of basic and non-basic variables.
	do {
		int coef = 0;
		for (int i = 0; i < var_num; ++i) {
			if (bitmask[i]) {
				for (int j = 0; j < res_num; j++) {
					matrix[j].push_back(res_matrix[j][i]);  // Build the matrix for the selected variables.
				}
				vars[coef] = i;
				coef++;
			}
		}

		// Step 4: Solve the linear system using the Gauss elimination method.
		solution = gauss(matrix, res_vector, res_num);

		// Step 5: Calculate the sum of the objective function coefficients for the current solution.
		double sum = 0;
		bool isAppr = false;

		// Check if any solution component is nonzero.
		for (int i = 0; i < res_num; i++) {
			if (solution[i] != 0) isAppr = true;
		}

		// Calculate the sum while considering infeasible solutions.
		for (int i = 0; i < res_num; i++) {
			if (solution[i] < 0) isAppr = false;  // Check for infeasibility.
			sum += solution[i] * fun_coeffs[vars[i]];  // Calculate the objective function value.
		}

		// Step 6: Check if the solution is both feasible and optimal.
		if (isAppr && sum < min_sum) {
			min_sum = sum;  // Update the minimum sum.
			best_solution = solution;  // Update the best solution found.
			vars_best = vars;  // Update the variables corresponding to the best solution.
		}
	} while (prev_permutation(bitmask.begin(), bitmask.end()));  // Generate all possible variable combinations.

	// Step 7: Create the final solution vector based on the best solution found.
	vector<double> sol(var_num);
	for (int i = 0; i < res_num; i++) {
		sol[vars_best[i]] = best_solution[i];
	}

	// Step 8: Return the final optimal solution.
	return sol;
}


vector<double> LPProblem::direct_by_dual(vector<double>& sol) {
	if (sol.size() != res_num)
		exit(-1);
	double eps = 0.00001;

	vector<int> nonzero;
	for (int i = 0; i < var_num; i++) {
		double sum = 0;
		for (int j = 0; j < res_num; j++) {
			sum += sol[j] * res_matrix[j][i];
		}
		if (abs(sum - fun_coeffs[i]) < eps) {
			nonzero.push_back(i);
		}
	}

	vector<vector<double>> matrix(nonzero.size());
	vector<double> vect(nonzero.size());
	int coef = 0;
	for (int i = 0; i < res_num; i++) {
		if (abs(sol[i]) > eps) {
			for (auto& cf : nonzero) {
				matrix[coef].push_back(res_matrix[i][cf]);
			}
			vect[coef] = res_vector[i];
			coef++;
		}
	}
	
	vector<double> gaus = gauss(matrix, vect, nonzero.size());
	vector<double> ans(var_num);
	coef = 0;
	for (auto& cf : nonzero) {
		ans[cf] = gaus[coef++];
	}
	
	return ans;
}

vector<double> LPProblem::gauss(vector<vector<double>>& A, vector<double>& Y, int n) {
	vector<double> x(n);  // The solution vector.
	vector<vector<double>> a = A;  // Copy of the coefficient matrix.
	vector<double> y = Y;  // Copy of the right-hand side vector.
	double max;
	int k, index;
	const double eps = 0.00001;  // A small value to check for zero coefficients.
	k = 0;

	// Step 1: Start the Gaussian elimination process.
	while (k < n) {
		max = abs(a[k][k]);  // Initialize max with the absolute value of the diagonal element.
		index = k;

		// Step 2: Find the maximum element in the current column below the diagonal.
		for (int i = k + 1; i < n; i++) {
			if (abs(a[i][k]) > max) {
				max = abs(a[i][k]);
				index = i;
			}
		}

		// Step 3: If the maximum element is smaller than a small threshold (eps), exit with the current solution.
		if (max < eps) {
			return x;
		}

		// Step 4: Swap the current row with the row containing the maximum element (partial pivoting).
		for (int j = 0; j < n; j++) {
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}

		// Swap the corresponding element in the right-hand side vector.
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;

		// Step 5: Perform row operations to create zeros below the diagonal element.
		for (int i = k; i < n; i++) {
			double temp = a[i][k];
			if (abs(temp) < eps) continue;  // Skip division by nearly zero elements.

			// Divide the row by the diagonal element to make it equal to 1.
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;

			if (i == k) continue;  // Skip the current row (no need to subtract itself).

			// Perform row operations to create zeros below the diagonal element in other rows.
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;  // Move to the next column.
	}

	// Step 6: Perform back-substitution to find the solution vector.
	for (k = n - 1; k >= 0; k--) {
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}

	// Step 7: Return the solution vector x.
	return x;
}


vector<vector<double>> LPProblem::inverse_matrix(vector<vector<double>>& A) {
	int n = A.size();  // Determine the size of the matrix (n x n).
	vector<vector<double>> E(n);  // Create an identity matrix (n x n).
	vector<vector<double>> a = A;  // Copy of the input matrix.

	// Step 1: Initialize the identity matrix E with 1's on the diagonal and 0's elsewhere.
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			E[i].push_back((i == j) ? 1 : 0);
		}
	}

	double max;
	int k, index;
	const double eps = 0.00001;  // A small value used for checking zero coefficients.
	k = 0;

	// Step 2: Start the Gauss-Jordan elimination process to reduce the original matrix A to the identity matrix.
	while (k < n) {
		// Step 3: Find the maximum element in the current column below the diagonal (partial pivoting).
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++) {
			if (abs(a[i][k]) > max) {
				max = abs(a[i][k]);
				index = i;
			}
		}

		// Step 4: If the maximum element is smaller than a small threshold (eps), exit with an error matrix.
		if (max < eps) {
			return { {0} };
		}

		// Step 5: Swap the current row and corresponding row in the identity matrix E.
		for (int j = 0; j < n; j++) {
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
			temp = E[k][j];
			E[k][j] = E[index][j];
			E[index][j] = temp;
		}

		// Step 6: Perform row operations to create zeros below the diagonal element in both matrices.
		for (int i = k; i < n; i++) {
			double temp = a[i][k];
			if (abs(temp) < eps) continue;  // Skip division by nearly zero elements.

			// Divide the current row by the diagonal element to make it equal to 1.
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] / temp;
				E[i][j] = E[i][j] / temp;
			}

			// Apply the same row operations to both matrices.
			if (i == k) continue;  // Skip the current row (no need to subtract itself).
			for (int j = 0; j < n; j++) {
				a[i][j] = a[i][j] - a[k][j];
				E[i][j] = E[i][j] - E[k][j];
			}
		}
		k++;  // Move to the next column.
	}

	// Step 7: Perform back-substitution on the identity matrix E to obtain the inverse matrix.
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = 0; k < n; k++) {
				E[j][k] -= E[i][k] * a[j][i];
			}
		}
	}

	// Step 8: Return the inverse matrix E.
	return E;
}
}