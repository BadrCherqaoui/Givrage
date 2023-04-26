#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>

#include "Mesh.h"
#include "operations.h"

using namespace std;


void lecture_obstacle(string fichier, vector<vector<double>> tab_noeud, vector<vector<int>> & tab_obstacle){

    int nb_triangles = 0; 
    int ignore_int;
    double ignore_double; 
    int nb_noeuds = 0; 
    int nb_cylinder = 0; 
    string ignore; 
    ifstream mon_flux; 
    mon_flux.open(fichier, ios::in); 
    
    mon_flux >> ignore >> ignore_int;
    mon_flux >> ignore >> nb_triangles;

    for (int i=0; i<nb_triangles; i++) {
        mon_flux >> ignore_int >> ignore_int >> ignore_int >> ignore_int >> ignore_int;
    }

    mon_flux >> ignore >> nb_noeuds;

    for (int i=0; i<nb_noeuds; i++) {
        mon_flux >> ignore_double >> ignore_double >> ignore_int;
    }

    getline(mon_flux, ignore); 
    getline(mon_flux, ignore);
    getline(mon_flux, ignore);
    mon_flux >> ignore >> nb_cylinder;
    tab_obstacle.resize(nb_cylinder, vector<int>(2)); 
    int noeud_min = -1;
    double s_min = 1000.;

    int i_n_min; 
    for (int i=0; i<nb_cylinder; i++) {
        mon_flux >> ignore_int >> tab_obstacle[i][0] >> tab_obstacle[i][1];
        if (s_min > abs(tab_noeud[tab_obstacle[i][0]][0]) + abs(tab_noeud[tab_obstacle[i][0]][1])){
            s_min = abs(tab_noeud[tab_obstacle[i][0]][0]) + abs(tab_noeud[tab_obstacle[i][0]][1]);
            noeud_min = tab_obstacle[i][0];
            i_n_min = i; 
        }
    }

    ofstream mon_flux_out; 
    mon_flux_out.open("Affichage/Obstacle.dat", ios::out); 

    
    for(int i = 0; i<nb_cylinder; i++){

        mon_flux_out <<  tab_noeud[tab_obstacle[i][0]][0] << " " <<  tab_noeud[tab_obstacle[i][0]][1] << endl; 
        mon_flux_out <<  tab_noeud[tab_obstacle[i][1]][0] << " " <<  tab_noeud[tab_obstacle[i][1]][1] << endl; 
    }

    mon_flux.close(); 

}


void lecture_maillage(string fichier, vector<vector<double>> & tab_noeud, vector<vector<int>> & tab_triangle, double & area_min, vector<vector<int>> & tab_voisin, vector<double> & tab_rho_a, vector<vector<double>> &tab_vitesse, double & v_max){

    int nb_triangles = 0; 
    int ignore_int; 
    int nb_noeuds = 0; 
    double area;
    string ignore; 
    ifstream mon_flux; 
    mon_flux.open(fichier, ios::in); 

    getline(mon_flux, ignore); 
    getline(mon_flux, ignore);
    getline(mon_flux, ignore); 
    getline(mon_flux, ignore); 
    
    mon_flux >> ignore >> nb_noeuds >> ignore;


    // Lecture des noeuds
    tab_noeud.resize(nb_noeuds, vector<double> (3));

    for (int i=0; i<nb_noeuds; i++) {
        mon_flux >> tab_noeud[i][0] >> tab_noeud[i][1] >> tab_noeud[i][2];
    }


    //Lecture triangles et calcul de l'aire min
    mon_flux >> ignore >> nb_triangles >> ignore_int;

    tab_triangle.resize(nb_triangles, vector<int> (3));

    area_min = 100000;
    for (int i=0; i<nb_triangles; i++) {
        mon_flux >> ignore_int >> tab_triangle[i][0] >> tab_triangle[i][1] >> tab_triangle[i][2];

        area = 0.5 * (tab_noeud[tab_triangle[i][0]][0]*(tab_noeud[tab_triangle[i][1]][1] - tab_noeud[tab_triangle[i][2]][1]) + tab_noeud[tab_triangle[i][1]][0]*(tab_noeud[tab_triangle[i][2]][1] - tab_noeud[tab_triangle[i][0]][1]) + tab_noeud[tab_triangle[i][2]][0]*(tab_noeud[tab_triangle[i][0]][1] - tab_noeud[tab_triangle[i][1]][1]) );
        
        area_min = min(area_min, area);
    }
  
    mon_flux >> ignore >> ignore_int; 

    //ignore 
    for (int i=0; i<nb_triangles; i++) {
        mon_flux >> ignore_int;
    }

    mon_flux >> ignore >> ignore_int;

    getline(mon_flux, ignore); 
    getline(mon_flux, ignore); 
    getline(mon_flux, ignore); 
    getline(mon_flux, ignore);


    //Lecture masses volumiques
    tab_rho_a.resize(nb_noeuds);

    for (int i=0; i<nb_noeuds; i++) {
        mon_flux >> tab_rho_a[i] ;
    }

    //lecture moments = rho*v
    getline(mon_flux, ignore);
    getline(mon_flux, ignore);
    tab_vitesse.resize(nb_noeuds, vector<double>(3));

    v_max = -1;

    for (int i=0; i<nb_noeuds; i++) {
        mon_flux >> tab_vitesse[i][0] >> tab_vitesse[i][1] >> tab_vitesse[i][2];

        for (int j=0; j<3; j++){
            tab_vitesse[i][j] = tab_vitesse[i][j]/tab_rho_a[i];
        }

        v_max = max(v_max, norm(tab_vitesse[i]));

    }

    mon_flux.close();


    // On cherche les voisins 
    vector<int> S(3,0); 
    tab_voisin.resize(nb_noeuds); 
    for(int i = 0; i<nb_triangles; i++){
        
        for(int j=0;j<3; j++){
            S[j] = tab_triangle[i][j]; 
            tab_voisin[S[j]].push_back(i); 
        }

    }

}


// Parcours tout un maillage afin de trouver le triangle ou se situe la particule d'eau
int where(vector<double> pos, vector<vector<int>> tab_triangle, vector<vector<double>> tab_noeud){

    int n1, n2, n3;  
    double sig1, sig2, sig3; 

    for(int i = 0; i<tab_triangle.size(); i++){

        n1 = tab_triangle[i][0]; 
        n2 = tab_triangle[i][1]; 
        n3 = tab_triangle[i][2]; 

        sig1 = (tab_noeud[n1][0] - pos[0])*(tab_noeud[n2][1] - pos[1]) - (tab_noeud[n1][1] - pos[1])*(tab_noeud[n2][0] - pos[0]); 
        sig2 = (tab_noeud[n2][0] - pos[0])*(tab_noeud[n3][1] - pos[1]) - (tab_noeud[n2][1] - pos[1])*(tab_noeud[n3][0] - pos[0]); 
        sig3 = (tab_noeud[n3][0] - pos[0])*(tab_noeud[n1][1] - pos[1]) - (tab_noeud[n3][1] - pos[1])*(tab_noeud[n1][0] - pos[0]); 

        if ( (sig1*sig2 >= 0) && (sig2*sig3 >= 0) ){
            return i; 
        }
    }

    cout << "Sortie de maillage" << endl;
    return -1; 

}

//Fonction permettant de chercher la particule uniquelent dans les triangles voisins du triangle où se trouvait la particule 
int where_voisin(int triangle_prec, std::vector<double> pos, std::vector<std::vector<int>> tab_triangle, std::vector<std::vector<double>> tab_noeud, std::vector<std::vector<int>> & tab_voisin, string fichier){

    vector<int> noeuds(3);

    //On obtient les noeuds du triangle ou se trouvait la particule 
    for(int i = 0; i<3; i++){
        noeuds[i] = tab_triangle[triangle_prec][i];
    }


    //Creation d'un tableau de numéros des voisins du triangles grace aux 3 noeuds du triangle
    vector<int> voisins(0); 
    for(int i = 0; i<3; i++){
        for(int j = 0; j<tab_voisin[noeuds[i]].size(); j++){
            if (!(count(voisins.begin(), voisins.end(), tab_voisin[noeuds[i]][j]))){
                voisins.push_back(tab_voisin[noeuds[i]][j]); 
            }
        }
    }

    // On cree donc les 3 noeuds d'un triangle voisins et on voit si la particule est dans un voisin 
    int n1, n2, n3; 
    double sig1, sig2, sig3; 

    for(int i = 0; i<voisins.size(); i++){

        n1 = tab_triangle[voisins[i]][0]; 
        n2 = tab_triangle[voisins[i]][1]; 
        n3 = tab_triangle[voisins[i]][2]; 
        
        sig1 = (tab_noeud[n1][0] - pos[0])*(tab_noeud[n2][1] - pos[1]) - (tab_noeud[n1][1] - pos[1])*(tab_noeud[n2][0] - pos[0]); 
        sig2 = (tab_noeud[n2][0] - pos[0])*(tab_noeud[n3][1] - pos[1]) - (tab_noeud[n2][1] - pos[1])*(tab_noeud[n3][0] - pos[0]); 
        sig3 = (tab_noeud[n3][0] - pos[0])*(tab_noeud[n1][1] - pos[1]) - (tab_noeud[n3][1] - pos[1])*(tab_noeud[n1][0] - pos[0]); 

        if ( (sig1*sig2 >= 0) && (sig2*sig3 >= 0) ){
            return voisins[i]; 
        }
    }

    // Sinon on est soit sorti du maillage soit on a sauté 2 mailles 
    if ((fichier == "mesh_tri_dim.su2") && (pos[0]>0.17)){
        cout << "Sortie du domaine" << endl;
        return -2; 
    }
    else if ((fichier == "meshtest.su2") && (pos[0]>0.6)){
        cout << "Sortie du domaine" << endl;
        return -2; 
    }
    else {
        cout << "Impact avec l'obstacle" << endl;
        return -1; 
    }

} 


