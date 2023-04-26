#ifndef OPERATIONS_H

#include <iostream>
#include <vector>
#include "operations.h"

std::vector<double> f(double t, std::vector<double> sol, double d, double rho, double rho_a, double mu_a, std::vector<double> ua, double Re, double Cd); 

void RK2(double t,std::vector<double> &sol,double dt, double d, double rho, double rho_a, double mu_a, std::vector<double> ua, double & Re); 

void RK4(double t,std::vector<double> &sol,double dt, double d, double rho, double rho_a, double mu_a, std::vector<double> ua, double & Re); 

double calcul_Re(std::vector<double> sol,double d, double rho, double rho_a, double mu_a, std::vector<double> ua);

#define OPERATIONS_H
#endif