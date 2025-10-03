#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <vector>
namespace la
{

void printMatrix(const std::vector<std::vector<double>>& matrix);
int gauss(std::vector<std::vector<double>> a, std::vector<double>& ans);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
}
#endif