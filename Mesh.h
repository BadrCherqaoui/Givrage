#include <vector>
#include <string>

void lecture_obstacle(std::string fichier, std::vector<std::vector<double>> tab_noeud, std::vector<std::vector<int>> & tab_obstacle); 

void lecture_maillage(std::string fichier, std::vector<std::vector<double>> & tab_noeud, std::vector<std::vector<int>> & tab_triangle, double & area_min, std::vector<std::vector<int>> & tab_voisin, std::vector<double>& tab_rho_a, std::vector<std::vector<double>>& tab_vitesse, double & v_max);

int where(std::vector<double> pos, std::vector<std::vector<int>> tab_triangle, std::vector<std::vector<double>> tab_noeud);

int where_voisin(int triangle_prec, std::vector<double> pos, std::vector<std::vector<int>> tab_triangle, std::vector<std::vector<double>> tab_noeud, std::vector<std::vector<int>> & tab_voisin, std::string fichier);

