#include <iostream>
#include <vector>
#include "LPProblem.h"

using namespace std;

int main() {
	setlocale(LC_ALL, "Russian");
	vector<double> f_coeffs = { 2, 1, 1 };

	vector<vector<double>> r_matrix = { {1,  2,  3},
										{2, 3,  1},
										{3,  1, 4} };
	vector<double> r_vector = { 1, 2, 3 };

	vector<int> sign = { 2, 0, 1 };

	vector<int> r_signs = { 1, 1, 0};

	cout << "Direct Task:\n\n";
	LPProblem* lp = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 0);
	lp->print();
	
	
	lp->to_canonical_form();
	cout << "Canonical form for Direct task:\n\n";
	lp->to_canonical_form();
	lp->print();
	vector<double> x = lp->simplex();
	cout << "Solution Simplex: " << x[0] << " " << x[1] << " " << x[2]-x[3] << "\n";
	cout << "Obective Function for direct task : " << lp->obj_fun(x) << "\n\n";
	x = lp->extreme_point_enum();
	cout << "Solution Extreme Points: " << x[0] << " " << x[1] << " " << x[2] - x[3] << "\n";
	cout << "Objection function for Direct task: " << lp->obj_fun(x) << "\n\n";

	cout << "Dual task:\n\n";
	LPProblem* lp2 = new LPProblem(f_coeffs, r_matrix, r_vector, sign, r_signs, 0);
	LPProblem* dual_lp = lp2->dual();
	dual_lp->print();
	cout << "Canonical form for Dual task:\n\n";
	dual_lp->to_canonical_form();
	dual_lp->print();
	x = dual_lp->simplex();
	
	cout << "Solution Simplex Method: " << x[0] << " " << x[1] << " " << x[2] - x[3] << "\n";
	cout << "Objection function for Dual task: " << dual_lp->obj_fun(x) << "\n\n";
	x = dual_lp->extreme_point_enum();
	cout << "Solution Extreme Points: " << x[0] << " " << x[1] << " " << x[2] - x[3] << "\n";
	cout << "Objection function for Dual task: " << dual_lp->obj_fun(x) << "\n\n";

}