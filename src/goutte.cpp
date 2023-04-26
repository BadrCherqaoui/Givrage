#ifndef GOUTTE_CPP

#include "goutte.h"
#include "operations.h"
#include "timescheme.h"
#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include "Mesh.h"




using namespace std;


//Générer des nombres random et uniforme 
std::default_random_engine generator(time(0)); 
std::uniform_real_distribution<double> distribution(-0.08, 0.08); 
std::uniform_real_distribution<double> distribution2(-40.0e-6, 40e-6); 
std::uniform_real_distribution<double> distribution3(-25.0, 25.0); 


goutte::goutte(int id): 
_id(id)
{

}


//==================== Initialisation de tous les parametres et ecrit dans un fichier les données d'entrée 
//===================== Ajout de condition aléatoire ======================================================
void goutte::initialisation(double d_min){

    string ignore;  
    
    bool rand = false; 
    double random_y = 0 ;
    double random_d = 0;
    double random_z = 0;

    double dim;

    if (rand){
        random_y = distribution(generator); 
        random_d= distribution2(generator); 
        random_z = distribution3(generator); 
    }
    
 
	std::ifstream monFlux("param_goutte.data", std::ios::in);


    if (monFlux){

        monFlux >> ignore; 
        monFlux >> dim;

        _pos_init.resize(dim); 
        _sol.resize(2*dim); 
        _ua.resize(dim); 

        monFlux >> ignore;


        if (dim == 2){
          monFlux >> _sol[0] >> _sol[1] >> _sol[2] >> _sol[3];

          _sol[1] += random_y;
          _sol[1] += _id * d_min * 6 ; // changer le d_min pour des maillages différents 
          _pos_init[0] = _sol[0];  
          _pos_init[1] = _sol[1];
        }
        else if(dim==3){
          monFlux >> _sol[0] >> _sol[1] >> _sol[2] >> _sol[3] >> _sol[4] >> _sol[5];
          
          _sol[1] += random_y;
          _sol[1] += _id * d_min;
          _pos_init[0] = _sol[0];  
          _pos_init[1] = _sol[1];
          _pos_init[2] = _sol[2];

        }

        else{
            cout << "Pas la bonne dimension"<< endl; 
            exit(1); 
        }

        if (dim == 2){
        }
        else if(dim==3){
        }
        else{
            cout << "Pas la bonne dimension"<< endl; 
            exit(1); 
        }
               

        monFlux >> ignore; 
        monFlux >> _rho ; 
        monFlux >> ignore; 
        monFlux >> _d; 
        monFlux >> ignore; 
        monFlux >> _t0; 
        monFlux >> ignore; 
        monFlux >> _dt; 
        monFlux >> ignore; 
        monFlux >> _tmax; 
        monFlux >> ignore; 
        monFlux >> _mu_a; 

        monFlux >> ignore; 
        monFlux >> _choix_interpolation;  

        monFlux >> ignore; 
        monFlux >> fichier_mesh; 

        monFlux >> ignore; 
        monFlux >> fichier_obstacle; 

        monFlux.close(); 


    }
    else{

        cout << "ERREUR: Impossible d'ouvrir le fichier en lecture" << endl;
        exit(1); 
    }


    _d += random_d;  
    if(dim==3){
      _sol[2] += random_z; 
    }

    

}


// ======================================================================================
// ==================== Différentes méthodes d'interpolation ============================
// ======================================================================================
void interpolation(vector<vector<double>> tab_noeud, vector<vector<int>> tab_triangle, vector<vector<double>> tab_vitesse, int i_triangle, vector<double> sol, vector<double> & ua, char choix_interpolation, double dmin, vector<double> tab_rho_a, double & rho_a){

    int dim = ua.size(); 
    vector<double> nodes(3, 0), pos_n1(3, 0), pos_n2(3, 0), pos_n3(3, 0);
    vector<double>  ua_n1(3,0),ua_n2(3,0),ua_n3(3,0); 
    vector<double> n1_M(dim, 0), n2_n1(dim, 0), n3_n1(dim, 0);
    vector<double>  rho_a_n(3,0);

    for (int i = 0; i<3; i++){
        nodes[i] =  tab_triangle[i_triangle][i];
        rho_a_n[i] = tab_rho_a[nodes[i]];
    }
            

    for (int i = 0; i<dim; i++){
        pos_n1[i] = tab_noeud[nodes[0]][i]; 
        pos_n2[i] = tab_noeud[nodes[1]][i]; 
        pos_n3[i] = tab_noeud[nodes[2]][i]; 

        ua_n1[i] = tab_vitesse[nodes[0]][i]; 
        ua_n2[i] = tab_vitesse[nodes[1]][i]; 
        ua_n3[i] = tab_vitesse[nodes[2]][i];

    }

    switch (choix_interpolation)
    {
    case 'm':

        for(int i = 0; i<dim; i++){
            ua[i] = (ua_n1[i]+ua_n2[i]+ua_n3[i])/3.;
        }

        rho_a = (rho_a_n[0] + rho_a_n[1] + rho_a_n[2])/3.;    

        break;

    case 'p':

        double delta1m, delta2m, delta3m, somme_delta; 

        for (int i = 0; i<dim; i++){
            delta1m += pow((pos_n1[i]-sol[i]),2); 
            delta2m += pow((pos_n2[i]-sol[i]),2); 
            delta3m += pow((pos_n3[i]-sol[i]),2); 
        }

        delta1m = sqrt(delta1m); 
        delta2m = sqrt(delta2m); 
        delta3m = sqrt(delta3m); 

        if (delta1m < dmin/100){
            for(int i = 0; i<dim; i++){
                ua[i] = ua_n1[i]; 
            }
            rho_a = rho_a_n[0];
        }
        else if (delta2m < dmin/100){
            for(int i = 0; i<dim; i++){
                ua[i] = ua_n2[i]; 
            }
            rho_a = rho_a_n[1];
        }
        else if (delta3m < dmin/100){
            for(int i = 0; i<dim; i++){
                ua[i] = ua_n3[i]; 
            }
            rho_a = rho_a_n[2];
        }
        else{

            somme_delta = 1./delta1m + 1./delta2m + 1./delta3m; 
        
            for(int i = 0; i<dim; i++){
                ua[i] = (1./delta1m) * ua_n1[i] / (somme_delta) + (1./delta2m) * ua_n2[i] / (somme_delta)+(1./delta3m) * ua_n3[i] / (somme_delta); 
            }

            rho_a = (1./delta1m) * rho_a_n[0] / (somme_delta) + (1./delta2m) * rho_a_n[1] / (somme_delta)+(1./delta3m) * rho_a_n[2] / (somme_delta); 
        
        }  

        break;

    case 'v' :

        double alpha, beta, det;

        for (int i=0; i<dim; i++) {
            n1_M[i] = sol[i] - pos_n1[i];
            n2_n1[i] = pos_n2[i] - pos_n1[i];
            n3_n1[i] = pos_n3[i] - pos_n1[i];
        }

        det = (n2_n1[0]*n3_n1[1] - n3_n1[0]*n2_n1[1]);

        alpha = (n3_n1[1]*n1_M[0] - n3_n1[0]*n1_M[1] ) / det;
        beta = (-n2_n1[1]*n1_M[0] + n2_n1[0]*n1_M[1] ) / det;

        for (int i=0; i<dim; i++) {
            ua[i] = (1-alpha-beta)*ua_n1[i] + alpha*ua_n2[i] + beta*ua_n3[i];
        }

        rho_a = (1-alpha-beta)*rho_a_n[0] + alpha*rho_a_n[1] + beta*rho_a_n[2];

        break;

    default:
        cout << "Choisissez une interpolation valide, nous avons utilisé la moyenne" << endl; 
        for(int i = 0; i<dim; i++){
            ua[i] = (ua_n1[i]+ua_n2[i]+ua_n3[i])/3.;  ; 
        }
        rho_a = (rho_a_n[0] + rho_a_n[1] + rho_a_n[2])/3.;    

        break;
    }
    

}



void find_edge(vector<vector<int>> tab_obstacle, int i_old, vector<vector<int>> tab_triangle, vector<vector<int>> & edges, vector<vector<double>> tab_noeud) {
    double e1, noeud_tr, noeud_tr2, noeud_tr3, e2;
    vector<int> edge(2,-1);
    int ind_1;
    int n_obst = tab_obstacle.size();

    noeud_tr = tab_triangle[i_old][0];
    noeud_tr2 = tab_triangle[i_old][1];
    noeud_tr3 = tab_triangle[i_old][2];

    //premier edge de tab_obstacle
    e1 = tab_obstacle[0][0];
    if ((e1 == noeud_tr) || (e1 == noeud_tr2) || (e1 == noeud_tr3)) {
        edge[0] = e1;
        ind_1 = 0;


        //conditions de bord
        e2 = tab_obstacle[0][1]; 
        if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
            edge[1] = e2;
        }
        else {
            // e2 = tab_obstacle[1][0];
            e2 = tab_obstacle[n_obst-1][0];
            if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
                edge[1] = e2;
            }
        }
    }
    

    for (int e=1; e<tab_obstacle.size()-1; e++) {
        e1 = tab_obstacle[e][0];
        if ((e1 == noeud_tr) || (e1 == noeud_tr2) || (e1 == noeud_tr3)) {
            edge[0] = e1;
            ind_1 = e;

            //conditions de bord
            e2 = tab_obstacle[e-1][0]; 
            if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
                edge[1] = e2;
            }
            else {
                e2 = tab_obstacle[e+1][0];
                if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
                    edge[1] = e2;
                }
            }
        }
    }

    //dernier edge de tab_obstacle
    e1 = tab_obstacle[n_obst-1][0];
    if ((e1 == noeud_tr) || (e1 == noeud_tr2) || (e1 == noeud_tr3)) {
        edge[0] = e1;
        ind_1 = tab_obstacle.size()-1;


        //conditions de bord
        e2 = tab_obstacle[tab_obstacle.size()-2][0]; 
        if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
            edge[1] = e2;
        }
        else {
            e2 = tab_obstacle[tab_obstacle.size()-1][1];
            if ((e2 == noeud_tr) || (e2 == noeud_tr2) || (e2 == noeud_tr3)) {
                edge[1] = e2;
            }
        }
    }




    //get the surroundings edges
    if (edge[1] == -1) {
        edges.push_back(tab_obstacle[ind_1]);

        if (ind_1 == n_obst-1) {
            edges.push_back(tab_obstacle[0]);
        }
        else {
            //discussion à avoir
            if (tab_noeud[edge[0]][1] < 0) {
                edges.push_back(tab_obstacle[ind_1 - 1]);
            }
            else {
                edges.push_back(tab_obstacle[ind_1 + 1]);
            }
        }
        
    } 
    else {
        if (ind_1 == 0) {
            edges.push_back(tab_obstacle[tab_obstacle.size()-1]);
        }
        else {
            edges.push_back(tab_obstacle[ind_1 - 1]);
        }

        edges.push_back(tab_obstacle[ind_1]);

        if (ind_1 == tab_obstacle.size()-1) {
            edges.push_back(tab_obstacle[0]);
        }
        else {
            edges.push_back(tab_obstacle[ind_1 + 1]);
        }
        
    }

}




void intersection(vector<vector<int>> edges, vector<double> old_pos, vector<double> pos, vector<vector<double>> tab_noeuds, vector<double> & impact){
    
    double a_pos, b_pos, a_e1, b_e1, a_e2, b_e2, a_e3, b_e3;
    vector <double> p1, p2, p3, p4;
    vector<vector<double>> pos_imp;

    // Droite position
    a_pos = (pos[1] - old_pos[1]) / (pos[0] - old_pos[0]);
    b_pos = (old_pos[1]*pos[0] - old_pos[0]*pos[1]) / (pos[0] - old_pos[0]);

    p1 = tab_noeuds[edges[0][0]];
    p2 = tab_noeuds[edges[1][0]];
    p3 = tab_noeuds[edges[1][1]];


    // Droite edge 1
    a_e1 = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    b_e1 = (p1[1]*p2[0] - p1[0]*p2[1]) / (p2[0] - p1[0]);

    // Droite edge 2
    a_e2 = (p3[1] - p2[1]) / (p3[0] - p2[0]);
    b_e2 = (p2[1]*p3[0] - p2[0]*p3[1]) / (p3[0] - p2[0]);

    // Droite edge 3 si il y a 
    if (edges.size() == 3)
    {
        p4 = tab_noeuds[edges[2][1]];
        a_e3 = (p4[1] - p3[1]) / (p4[0] - p3[0]);
        b_e3 = (p3[1]*p4[0] - p3[0]*p4[1]) / (p4[0] - p3[0]);
    }

    pos_imp.resize(3, {-1,-1});
    // premiere intersection x et y
    pos_imp[0][0] = (b_e1 - b_pos) / (a_pos - a_e1);
    pos_imp[0][1] = a_pos * pos_imp[0][0] + b_pos;
    
    // deuxieme intersection x et y
    pos_imp[1][0] = (b_e2 - b_pos) / (a_pos - a_e2);
    pos_imp[1][1] = a_pos * pos_imp[1][0] + b_pos;

    if (edges.size() == 3)
    {
        pos_imp[2][0] = (b_e3 - b_pos) / (a_pos - a_e3);
        pos_imp[2][1] = a_pos * pos_imp[2][0] + b_pos;
    }

    if ( (pos_imp[0][0] <= max(p2[0], p1[0])) && (pos_imp[0][0] >= min(p2[0], p1[0])) && (pos_imp[0][1] <= max(p2[1], p1[1])) && (pos_imp[0][1] >= min(p2[1], p1[1])) ) {
        impact[0] = pos_imp[0][0];
        impact[1] = pos_imp[0][1];
    } 

    else if ( (pos_imp[1][0] <= max(p2[0], p3[0])) && (pos_imp[1][0] >= min(p2[0], p3[0])) && (pos_imp[1][1] <= max(p2[1], p3[1])) && (pos_imp[1][1] >= min(p2[1], p3[1])) ) {
        impact[0] = pos_imp[1][0];
        impact[1] = pos_imp[1][1];
    } 

    else if (edges.size() == 3) {
        if ( (pos_imp[2][0] <= max(p4[0], p3[0])) && (pos_imp[2][0] >= min(p4[0], p3[0])) && (pos_imp[2][1] <= max(p4[1], p3[1])) && (pos_imp[2][1] >= min(p4[1], p3[1])) ) {
        impact[0] = pos_imp[2][0];
        impact[1] = pos_imp[2][1];
        }
    }

    else {
        cout << "VOIR INTERSECTION D'URGENCE" << endl;
    }

}





// ================================= Fonction avancement total d'une particule =================================
void goutte::advance(std::vector<std::vector<int>> tab_obstacle, std::vector<std::vector<double>> tab_noeud, std::vector<std::vector<int>> tab_triangle, std::vector<double> tab_rho_a, std::vector<std::vector<double>> tab_vitesse, std::vector<std::vector<int>> tab_voisin, double area_min, double vitesse_max){

    int dim=_sol.size()/2;    
    double Re; 

    int every_points = 10;

    vector<double> old_pos(dim, 0); 
    _pos_imp.resize(dim); 

    vector<int> nodes(3,0); 
    vector<double> ua_0(dim,0);
    vector<double> ua_n1(dim,0), ua_n2(dim,0), ua_n3(dim,0); 
    vector<double> pos_n1(dim,0), pos_n2(dim,0), pos_n3(dim,0); 


    // =============================== Intégration des éléments du Maillage ===============================
    int i_triangle, i_old; 
    double d_min, dt_max; 


    d_min = sqrt(area_min);

    dt_max = d_min/vitesse_max; 
    if (dt_max < _dt){
        _dt = 0.5*min(dt_max,1e-4); 
        cout << "------------------- Changement de pas de temps, nous avons pris dt = " << _dt << "-------------------"<< endl; 
    }


    // ======================================= Fichier résultat  =======================================
    ofstream mon_flux; 


    if (dim==2){
        mon_flux << "# t     x        y       vx      vy       Re" << endl; 
    }
    else if (dim==3){
        mon_flux << "# t     x        y       z      vx       vy      vz      Re" << endl; 
    }
    else{
        cout << "Pas la bonne dimension" << endl; 
        exit(1); 
    }

    // ============================================ Boucle temporelle =================================


    // ==================== On trouve où est la particule au début en parcourant tout le maillage  =====
    i_triangle = where(_sol, tab_triangle, tab_noeud);
    _t = _t0; 

    int n_iter = floor(_tmax/_dt) + 1; 
    int n_impact = 0; 
    
    vector<vector<double>> save_sol; 
    save_sol.resize(n_iter/every_points+1, vector<double> (2*dim+2)); 
    for(int k = 0; k<n_iter; k++){


        if (i_triangle >= 0){

            interpolation(tab_noeud, tab_triangle, tab_vitesse, i_triangle, _sol, _ua, _choix_interpolation, d_min, tab_rho_a, _rho_a); 
            
            for (int i = 0; i <dim; i++){
                old_pos[i] = _sol[i]; 
            }

            RK4(_t, _sol, _dt, _d, _rho, _rho_a, _mu_a, _ua, Re);

            if(k%every_points == 0){
                
                save_sol[k/every_points][0] = _t+_dt;
                for (int i = 0; i<2*dim; i++){
                    save_sol[k/every_points][i+1] = _sol[i]; 
                }
                save_sol[k/every_points][2*dim+1] =  Re; 

            }

            //On regarde ou est la particule au temps suivant grace au numéros du triangle d'avant 
            i_old = i_triangle;
            i_triangle = where_voisin(i_triangle, _sol, tab_triangle, tab_noeud, tab_voisin, fichier_obstacle); 
            _t += _dt;  

        }
        else if (i_triangle == -2){
            k = n_iter;  
        }
        else{
            n_impact = k;
            k = n_iter;  
        }

        n_impact = k; 
    }
    
    mon_flux.open("../Affichage/Resultats_"+to_string(_id)+ "_"+ _choix_interpolation + ".dat", ios::out);
    for(int i = 0; i<(n_impact/every_points); i++){
        for(int j = 0;j<2*dim+2; j++){
            mon_flux << save_sol[i][j] << " "; 
        }
        mon_flux << " " << endl; 
    }
    mon_flux.close(); 


    if (_choix_interpolation == 'm'){
        cout << "------------------- Utilisation de la méthode d'interpolation moyenne -------------------" << endl; 
    }
    else if(_choix_interpolation == 'v'){
        cout << "------------------- Utilisation de la méthode d'interpolation vecteur -------------------" << endl; 
    }
    else if(_choix_interpolation == 'p'){
        cout << "------------------- Utilisation de la méthode d'interpolation poids -------------------" << endl; 

    }
    else{
        cout << "------------------- Utilisation de la méthode d'interpolation moyenne -------------------" << endl; 
    }

    //====================== Si Impact ======================
    if (i_triangle == -1){

        vector<double> pos(2,0);
        vector<vector<int>> edges;
        vector<double> impact;
        impact.resize(2, 0);

        pos[0] = _sol[0];
        pos[1] = _sol[1]; 
     
        find_edge(tab_obstacle, i_old, tab_triangle, edges, tab_noeud);
        intersection(edges, old_pos, pos, tab_noeud, impact);

        cout << "Impact numérique :   " << impact[0] << " " << impact[1] << endl;

        std::ofstream mon_flux2;
        mon_flux2.open("../Affichage/Impact_"+std::to_string(_id)+".dat", std::ios::out);
        for(int i = 0; i<dim; i++){
            _pos_imp[i] = impact[i];
            mon_flux2 << _pos_imp[i] << " ";  
        }
        mon_flux2.close();
        
    }

}



#define GOUTTE_CPP
#endif
