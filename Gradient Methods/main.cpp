#include <iostream>
#include <cmath>
#include<vector>
using namespace std;


// Objective function: f(x, y, z) = (x - 2)^2 + (y - 3)^2 + (z - 4)^2
double objectiveFunction(double x, double y, double z) {
    return (x - 2) * (x - 2) + (y - 3) * (y - 3) + (z - 4) * (z - 4);
}

// Gradient of the objective function
void gradient(double x, double y, double z, double& grad_x, double& grad_y, double& grad_z) {
    grad_x = 2 * (x - 2);
    grad_y = 2 * (y - 3);
    grad_z = 2 * (z - 4);
}

void hess(double x, double y, double z, vector<vector<double>>& hess) {
    // Compute the second derivatives
    hess[0][0] = 2; // ∂²f/∂x²
    hess[1][1] = 2; // ∂²f/∂y²
    hess[2][2] = 2; // ∂²f/∂z²

    // Cross-derivatives (off-diagonal elements)
    hess[0][1] = 0; // ∂²f/∂x∂y
    hess[1][0] = 0; // ∂²f/∂y∂x
    hess[0][2] = 0; // ∂²f/∂x∂z
    hess[2][0] = 0; // ∂²f/∂z∂x
    hess[1][2] = 0; // ∂²f/∂y∂z
    hess[2][1] = 0; // ∂²f/∂z∂y
}


vector<vector<double>> inverse_matrix(vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> E(n);
    vector<vector<double>> a = A;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            E[i].push_back(((i == j) ? 1 : 0));
        }
    }
    double max;
    int k, index;
    const double eps = 0.000001;
    k = 0;
    while (k < n) {
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > max) {
                max = abs(a[i][k]);
                index = i;
            }
        }
        if (max < eps) {
            return { {0} };
        }
        for (int j = 0; j < n; j++) {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
            temp = E[k][j];
            E[k][j] = E[index][j];
            E[index][j] = temp;
        }
        for (int i = k; i < n; i++) {
            double temp = a[i][k];
            if (abs(temp) < eps) continue;
            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] / temp;
                E[i][j] = E[i][j] / temp;
            }
            if (i == k)  continue;
            for (int j = 0; j < n; j++) {
                a[i][j] = a[i][j] - a[k][j];
                E[i][j] = E[i][j] - E[k][j];
            }
        }
        k++;
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            for (int k = 0; k < n; k++) {
                E[j][k] -= E[i][k] * a[j][i];
            }
        }
    }
    return E;
}


// Gradient Descent with Step Size Reduction
void gradientDescentWithStepSizeReduction(double epsilon, double alpha0, double lambda, double delta) {
    double x = 0.0, y = 0.0, z = 0.0;
    int functionCalls = 0;
    double alpha = alpha0;

    while (true) {
        double grad_x, grad_y, grad_z;
        gradient(x, y, z, grad_x, grad_y, grad_z);
        functionCalls++;

        double f_current = objectiveFunction(x, y, z);
        double x_next = x - alpha * grad_x;
        double y_next = y - alpha * grad_y;
        double z_next = z - alpha * grad_z;
        double f_next = objectiveFunction(x_next, y_next, z_next);

        while (f_next - f_current > -alpha * delta * (grad_x * grad_x + grad_y * grad_y + grad_z * grad_z)) {
            alpha *= lambda;
            x_next = x - alpha * grad_x;
            y_next = y - alpha * grad_y;
            z_next = z - alpha * grad_z;
            f_next = objectiveFunction(x_next, y_next, z_next);
        }

        x = x_next;
        y = y_next;
        z = z_next;

        if (sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z) < epsilon) {
            break;
        }
    }

    std::cout << "Gradient Descent with Step Size Reduction:\n";
    std::cout << "Minimum point: (" << x << ", " << y << ", " << z << ")\n";
    std::cout << "Minimum value: " << objectiveFunction(x, y, z) << "\n";
    std::cout << "Function calls: " << functionCalls << "\n";
}

// Newton's Second-Order Gradient Descent

void newtonsSecondOrderGradientDescent(double epsilon, vector<vector<double>>& hessian) {
    double x = 0.0, y = 0.0, z = 0.0;
    int functionCalls = 0;

    while (true) {
        double grad_x, grad_y, grad_z;
        gradient(x, y, z, grad_x, grad_y, grad_z);
        hess(x, y, z, hessian);
        functionCalls++;

        double gradNorm = grad_x * grad_x + grad_y * grad_y + grad_z * grad_z;
        if (gradNorm < epsilon) {
            break;
        }

        vector<vector<double>> inv_hess = inverse_matrix(hessian); // Invert the Hessian

        // Compute the step using the inverted Hessian and the gradient
        double step_x = -(inv_hess[0][0] * grad_x + inv_hess[0][1] * grad_y + inv_hess[0][2] * grad_z);
        double step_y = -(inv_hess[1][0] * grad_x + inv_hess[1][1] * grad_y + inv_hess[1][2] * grad_z);
        double step_z = -(inv_hess[2][0] * grad_x + inv_hess[2][1] * grad_y + inv_hess[2][2] * grad_z);

        x += step_x;
        y += step_y;
        z += step_z;
    }

    std::cout << "Newton's Second-Order Gradient Descent:\n";
    std::cout << "Minimum point: (" << x << ", " << y << ", " << z << ")\n";
    std::cout << "Minimum value: " << objectiveFunction(x, y, z) << "\n";
    std::cout << "Function calls: " << functionCalls << "\n";
}


int main() {
    double epsilon = 1e-6;
    double alpha0 = 0.1;
    double lambda = 0.5;
    double delta = 0.1;

    vector<vector<double>> hessian(3, vector<double>(3));
    hess(0.0, 0.0, 0.0, hessian); // Calculate the Hessian matrix

    gradientDescentWithStepSizeReduction(epsilon, alpha0, lambda, delta);
    newtonsSecondOrderGradientDescent(epsilon, hessian);

    return 0;
}

