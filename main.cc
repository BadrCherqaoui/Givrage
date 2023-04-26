#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>
#include "captation.h"
#include "Mesh.h"
#include "goutte.h"
#include "mpi.h"


using namespace std;

// ======================================================================================
// ====================  Fonction de répartition de charge ==============================
// ======================================================================================
void charge(int me, int n, int np, int *iBeg, int *iEnd)
{
  int r = n%np;
  if (me<r)
    {
      *iBeg = me*(n/np) + me; 
      *iEnd = (me+1)*(n/np +1) - 1;
    }
  else
    {
      *iBeg = r + me * (n/np);
      *iEnd = *iBeg + n/np - 1;
    }
}



int main(int argc, char** argv){

    // =======================================================================================================
    // ===========================================  Creation du dossier Affichage ===========================
    // =======================================================================================================
    mkdir("Affichage", 0777);
    mkdir("Resultats", 0777);
    remove("Resultats/beta_en_fonction_de_s");


    // =========================================== Calcul de temps =========================================== 
    auto begin = std::chrono::high_resolution_clock::now();
    


    // =========================================== Choix nombre de gouttes ===================================
    int n_droplets(2);
    goutte* first_droplet;
    first_droplet = new goutte(-1);



    // ===========================================   Lecture du maillage  =====================================
    first_droplet -> initialisation(1.); //Lecture pour trouver quel fichier utiliser
    
    vector<vector<double>>  tab_vitesse, tab_noeud; 
    vector<double> tab_rho_a;
    vector<vector<int>> tab_voisin,tab_obstacle, tab_triangle; 
    double area_min, d_min, vitesse_max, dt_max;
    string fichier_mesh = first_droplet -> get_fichier(); 
    string fichier_obstacle = first_droplet -> get_fichier_obstacle(); 

    lecture_maillage(fichier_mesh, tab_noeud, tab_triangle, area_min, tab_voisin, tab_rho_a, tab_vitesse, vitesse_max); 
    lecture_obstacle(fichier_obstacle, tab_noeud, tab_obstacle); 



    //  =========================================== Tab goutte pour la BETA courbe  ===================================
    vector<vector<double>> tab_goutte; 
    tab_goutte.resize(n_droplets, vector<double> (4)); 

    int rank, np, iBeg, iEnd, s_loc; 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    charge(rank, n_droplets, np, &iBeg, &iEnd);
    s_loc = iEnd - iBeg + 1;

    vector<vector<double>> tab_goutte_loc; 
    tab_goutte_loc.resize(s_loc, vector<double> (4)); 

    
    d_min = sqrt(area_min)/2;
    
    for(int i = iBeg; i<iEnd+1 ;i++){
        cout << "------------------- Goutte " << i+1 << "/" << n_droplets << " -------------------" << endl;
        first_droplet = new goutte(i);
        first_droplet -> initialisation(d_min); 
        first_droplet -> advance(tab_obstacle, tab_noeud, tab_triangle, tab_rho_a, tab_vitesse, tab_voisin, area_min, vitesse_max);  
        tab_goutte_loc[i-iBeg][0] = first_droplet -> get_pos_init()[0];
        tab_goutte_loc[i-iBeg][1] = first_droplet -> get_pos_init()[1];
        tab_goutte_loc[i-iBeg][2] = first_droplet -> get_pos_imp()[0];
        tab_goutte_loc[i-iBeg][3] = first_droplet -> get_pos_imp()[1];
    }

    if (rank == 0) {

        int iBeg2, iEnd2, s_loc2;

        for (int i=0; i<s_loc; i++) {
            for(int j=0; j<4; j++){
                tab_goutte[i][j] = tab_goutte_loc[i][j];
            }
        }
        for (int r=1; r<np; r++) {
            charge(r, n_droplets, np, &iBeg2, &iEnd2);
            s_loc2 = iEnd2 - iBeg2 + 1;
            for(int i = iBeg2; i<iEnd2+1; i++){
                MPI_Recv(&tab_goutte[i][0], 4, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

    }
    else {

        for(int i = iBeg; i<iEnd+1; i++){
            MPI_Send(&tab_goutte_loc[i-iBeg][0],4, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

    }

    if (rank == 0) {

        //  =========================================== Calcul de BETA  ==============================================
        recherche(tab_goutte, tab_obstacle, tab_noeud); 


        //  =========================================== Fin de calcul du temps  ======================================
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);

    }

    MPI_Finalize(); 

    return 0;
}