#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

bool in_cv(std::vector<double> noeud, std::vector<double> goutte1, std::vector<double> goutte2); 
void recherche(std::vector<std::vector<double>> tab_goutte,std::vector<std::vector<int>> tab_obstacle, std::vector<std::vector<double>> tab_noeud); 