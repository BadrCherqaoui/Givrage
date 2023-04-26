#ifndef OPERATIONS_H

#include <iostream>
#include <vector>



bool in_obstacle(std::vector<double> sol, double R, std::vector<double> C0, double h);

std::vector<double> addition(std::vector<double> u, std::vector<double> v);

std::vector<double> prod_const_vector(double alpha, std::vector<double> u);

std::vector<double> moins(std::vector<double> u, std::vector<double> v);

double norm(std::vector<double> u);

double ps(std::vector<double> u, std::vector<double> v);




#define OPERATIONS_H
#endif
