#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

#include "captation.h"


using namespace std; 



// ====================================================================================================
// ====================  Fonction volume de controle pour impact numérque =============================
// ====================================================================================================
bool in_cv(vector<double> noeud, vector<double> goutte1, vector<double> goutte2){

    double xnoeud = noeud[0]; 
    double ynoeud = noeud[1]; 

    double xmin = min(goutte1[0],goutte2[0]); 
    double xmax = max(goutte1[0],goutte2[0]); 
    double ymax = max(goutte1[1],goutte2[1]); 
    double ymin = min(goutte1[1],goutte2[1]); 


    return ( ((ymin<ynoeud) && (ynoeud<ymax)) && ((xnoeud>xmin) && (xnoeud<xmax)) ) ; 

}



// ======================================================================================
// ====================  Recherche des impacts de gouttes pour le calcul de Beta ========
// ======================================================================================
void recherche(vector<vector<double>> tab_goutte,vector<vector<int>> tab_obstacle, vector<vector<double>> tab_noeud){

    int n_obst = tab_obstacle.size();
    int n_goutte = tab_goutte.size();  
    int noeud_etudie, noeud_present, noeud_suivant; 
    vector<double> noeud(2,0); 
    vector<double> goutte1(2,0), goutte2(2,0); 

    double dy, ds, beta, s=0; 

    ofstream mon_flux, mon_flux_entoure; 
    mon_flux.open("Resultats/beta_en_fonction_de_s", ios::app); 
   

    for (int j=0; j<n_goutte-1; j++){

        goutte1[0] = tab_goutte[j][2]; 
        goutte1[1] = tab_goutte[j][3]; 
        goutte2[0] = tab_goutte[j+1][2]; 
        goutte2[1] = tab_goutte[j+1][3]; 

        for (int i=0; i<n_obst; i++){
            noeud_etudie = tab_obstacle[i][0];
            noeud[0] = tab_noeud[noeud_etudie][0]; 
            noeud[1] = tab_noeud[noeud_etudie][1]; 

            if (in_cv(noeud, goutte1, goutte2)){
                dy = abs(tab_goutte[j][1]-tab_goutte[j+1][1]); 
                ds = sqrt(pow(goutte1[0]-goutte2[0],2)+pow(goutte1[1]-goutte2[1],2)); 
                cout << "dy = " << dy << "   ds = " << ds << endl;
                beta = dy/ds; 
                if (noeud[1]<0){
                    mon_flux << -sqrt(pow(noeud[0],2)+pow(noeud[1],2)) << " " << beta << endl; 
                }
                else{
                    mon_flux << sqrt(pow(noeud[0],2)+pow(noeud[1],2)) << " " << beta << endl; 
                }

                i = n_obst; 
            }

        }
    }

    mon_flux.close(); 

}
