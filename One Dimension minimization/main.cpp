#include <iostream>
#include <cmath>
using namespace std;

// Define the function to minimize
double function(double x) {
	return -x * x - x - sin(x);
}

// Define the function to minimize for the Fibonacci search
double f(double x) {
	return x * x + x + sin(x);
}

// Implement the Golden Ratio search algorithm
void GoldenRatio(double left, double right, double eps) {
	const double phi = (3 - sqrt(5)) / 2;
	double y, z, Fy, Fz;
	int funcCalls = 0;

	// Initialize the initial points and function values
	y = left + phi * (right - left);
	Fy = function(y);
	funcCalls++;

	z = right + left - y;
	Fz = function(z);
	funcCalls++;

	while (true) {
		if (Fy > Fz) {
			// Update the right endpoint and check for termination
			right = z;
			if (right - left < eps) break;

			// Update the other variables
			z = y;
			Fz = Fy;
			y = right + left - y;
			Fy = function(y);
			funcCalls++;
		}
		else {
			// Update the left endpoint and check for termination
			left = y;
			if (right - left < eps) break;

			// Update the other variables
			y = z;
			Fy = Fz;
			z = left + right - z;
			Fz = function(z);
			funcCalls++;
		}
	}

	// Output the result and number of function calls
	cout << "Golden Ratio Method:\nExtrem Point: " << (right + left) / 2. << " +- " << eps / 2 << endl;
	cout << "Function calls: " << funcCalls << "\n\n";
}

// Function to perform Fibonacci search for finding the minimum of a function.
double fibonacci_search(double (*f)(double), double a, double b, int n, double epsilon) {
	// Define constants phi and s based on the Fibonacci sequence
	double s = (1 - sqrt(5)) / (1 + sqrt(5));
	double phi = (1 + sqrt(5)) / 2.0;

	// Initialize variables for the Fibonacci search
	double rho = 1 / (phi * (1 - pow(s, n + 1)) / (1 - pow(s, n)));
	double d = rho * b + (1 - rho) * a;
	double yd = f(d);  // Evaluate the function at point d

	int funcCalls = 0;  // Keep track of the number of function evaluations

	// Perform the Fibonacci search for 'n' iterations
	for (int i = 1; i <= n - 1; ++i) {
		double c;
		double yc;
		if (i == n - 1) {
			// Calculate the final point 'c' using epsilon
			c = epsilon * a + (1 - epsilon) * d;
		}
		else {
			// Update the rho value based on the Fibonacci sequence
			rho = 1 / (phi * (1 - pow(s, n - i + 1)) / (1 - pow(s, n - i)));
			// Calculate the next point 'c'
			c = rho * a + (1 - rho) * b;
		}
		// Evaluate the function at point 'c'
		yc = f(c);
		funcCalls++;

		// Update the values a, b, d, and yd based on the function values
		if (yc < yd) {
			b = d;
			d = c;
			yd = yc;
		}
		else {
			a = b;
			b = c;
		}
	}

	// Output the result, including the estimated extremum point and function call count
	std::cout << "Fibonacci Method:\n Extrem Point: " << d << " +- " << epsilon / 2 << " \n Function calls: " << funcCalls << std::endl;
	return d;  // Return the estimated extremum point
}

// Function to find the minimum of a function within an interval using a simple method.
void findMinimum(double a, double b, double epsilon, int n) {
	// Calculate the initial step size 'h' based on the number of subdivisions 'n'
	double h = (b - a) / n;
	int funcCalls = 0;  // Keep track of the number of function evaluations

	// Perform the search for minimum
	while (h > epsilon) {
		double x1 = a;
		double x2 = a + h;

		// Keep moving along the interval until the function values decrease
		while (f(x1) > f(x2)) {
			x1 = x2;
			x2 = x1 + h;
			funcCalls++;
			if (x1 == b) {
				break;  // Exit if we reach the end of the interval
			}
		}

		// Update the interval 'a' and 'b' and recalculate the step size 'h'
		a = x1 - h;
		b = x2;
		h = (b - a) / n;
	}

	// Output the result, including the estimated extremum point and function call count
	std::cout << "Find Minimum:\nExtrem Point: " << (a + b) / 2. << " +- " << epsilon / 2 << std::endl;
	std::cout << "Function calls: " << funcCalls << "\n\n";
}

int main() {
	cout.setf(std::ios_base::fixed);
	cout.precision(7);

	const double LEFT = -2;
	const double RIGHT = 0;
	const double epsilon = 0.000001;
	int n = 10; // Initial number of subdivisions Начальное количество разбиений

	// Call the findMinimum, Fibonacci search, and Golden Ratio methods
	findMinimum(LEFT, RIGHT, epsilon, n);
	fibonacci_search(f, LEFT, RIGHT, n, epsilon);
	//GoldenRatio(LEFT, RIGHT, epsilon);  // Golden Ratio method is commented out

	return 0;
}
