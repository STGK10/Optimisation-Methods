#include <iostream>
#include <cmath>

// Определение целевой функции f(x, y, z)
double targetFunction(double x, double y, double z) {
    return x * x + 2 * y * y + 3 * z * z;
}

// Функция для расчета штрафной функции P(x, y, z, mu)
double penaltyFunction(double x, double y, double z, double mu1, double mu2, double mu3) {
    return targetFunction(x, y, z) + mu1 * std::pow(x + y - 5, 2) + mu2 * std::pow(std::max(0.0, 10 - x * y * z), 2) + mu3 * std::pow(std::max(0.0, z * z - 4 * x * y - 15), 2);
}

// Градиент штрафной функции P(x, y, z, mu) по переменным (x, y, z)
void gradient(double x, double y, double z, double mu1, double mu2, double mu3, double& gradX, double& gradY, double& gradZ) {
    gradX = 2 * x + 2 * mu1 * (x + y - 5) - 2 * mu2 * (z >= 10 ? 0 : -2 * y * z);
    gradY = 4 * y + 2 * mu1 * (x + y - 5) - 2 * mu2 * (z >= 10 ? 0 : -2 * x * z);
    gradZ = 6 * z + 2 * mu2 * (z >= 10 ? 0 : -x * y) - 2 * mu3 * (z * z - 4 * x * y - 15 >= 0 ? 0 : 2 * z);
}

// Метод градиентного спуска
void gradientDescent(double& x, double& y, double& z, double mu1, double mu2, double mu3, double learningRate, int maxIterations) {
    for (int i = 0; i < maxIterations; ++i) {
        double gradX, gradY, gradZ;
        gradient(x, y, z, mu1, mu2, mu3, gradX, gradY, gradZ);

        x -= learningRate * gradX;
        y -= learningRate * gradY;
        z -= learningRate * gradZ;
    }
}

// Метод золотого сечения для оптимизации коэффициентов штрафа
double goldenSectionSearch(double x, double y, double z, double mu, double epsilon) {
    double a = 0.0;
    double b = 1.0;
    double phi = (1.0 + std::sqrt(5.0)) / 2.0;

    double x1 = a + (1.0 - 1.0 / phi) * (b - a);
    double x2 = a + 1.0 / phi * (b - a);

    double f1 = penaltyFunction(x, y, z, x1, x1, x1);
    double f2 = penaltyFunction(x, y, z, x2, x2, x2);

    while (std::abs(b - a) > epsilon) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            x1 = a + (1.0 - 1.0 / phi) * (b - a);
            f2 = f1;
            f1 = penaltyFunction(x, y, z, x1, x1, x1);
        }
        else {
            a = x1;
            x1 = x2;
            x2 = a + 1.0 / phi * (b - a);
            f1 = f2;
            f2 = penaltyFunction(x, y, z, x2, x2, x2);
        }
    }

    return (a + b) / 2.0;
}

int main() {
    double x = 0.0, y = 0.0, z = 0.0; // Начальные значения переменных
    double mu1 = 1.0, mu2 = 1.0, mu3 = 1.0; // Начальные значения коэффициентов штрафа
    double learningRate = 0.1; // Скорость обучения
    int maxIterations = 100; // Максимальное количество итераций
    double epsilon = 1e-5; // Точность для метода золотого сечения

    gradientDescent(x, y, z, mu1, mu2, mu3, learningRate, maxIterations);

    // Оптимизация коэффициентов штрафа с использованием метода золотого сечения
    mu1 = goldenSectionSearch(x, y, z, mu1, epsilon);
    mu2 = goldenSectionSearch(x, y, z, mu2, epsilon);
    mu3 = goldenSectionSearch(x, y, z, mu3, epsilon);

    std::cout << "Optimal solution: x = " << x << ", y = " << y << ", z = " << z << std::endl;
    std::cout << "Optimal value: " << targetFunction(x, y, z) << std::endl;
    std::cout << "Optimal mu1: " << mu1 << ", mu2: " << mu2 << ", mu3: " << mu3 << std::endl;

    return 0;
}
