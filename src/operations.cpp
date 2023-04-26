#ifndef OPERATIONS_CPP

#include "timescheme.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


//==========================-  Calcul de norme ==========================
double norm(vector<double> u){
    double norm=0.;

    for(int i=0;i<u.size();i++)
    {
        norm+=u[i]*u[i];
    }
    norm=sqrt(norm);

    return norm;
}


//========================== Calcul de produit scalaire ==========================

double ps(vector<double> u, vector<double> v) {
    double prod_scalaire(0.);

    if (u.size() != v.size()) {
        cout << "VARIABLES DU PROD SCAL N'ONT PAS LES MEMES DIMENSIONS" << endl;
        exit(1);
    }

    for (int i=0; i<u.size(); i++) {
        prod_scalaire += u[i] * v[i]; 
    } 

    return prod_scalaire;
}


// ==========================  Si le point est dans le cercle ou le cylindre ========================
bool in_obstacle(vector<double> sol, double R, vector<double> C0, double h){
    int dim=sol.size()/2;

    vector<double> pos_xy(2,0);
    for(int i=0;i<2;i++) pos_xy[i]=sol[i];

    if (dim == 2)
    {
        return ((norm(moins(pos_xy,C0))<=R));
    }
    else if(dim==3){
        return ((norm(moins(pos_xy,C0))<=R) && (sol[2] <= h) && (sol[2] >= 0));
    }
    else{
        cout << "Pas la bonne dimension"<< endl; 
        exit(1); 
    }
}

// ====================================== Fonction opération vecteur ======================================

vector<double> addition(vector<double> u, vector<double> v){

    int n = v.size();
    int m = u.size();
    vector<double> res;
    res.resize(n);

    if (m==n){
        for(int i = 0; i<n;i ++){
            res[i] = u[i] + v[i];
        }
        return res;
    }
    else{
        cout << "les vecteurs u et v n'ont pas la mêmes dimensions pour le +"<< endl;
        return res;
    }

}

vector<double> prod_const_vector(double alpha, vector<double> u){

    vector<double> v;
    int n = u.size();
    v.resize(n) ;

    for(int i = 0; i<n; i++){
        v[i] = alpha*u[i];
    }

    return v;
}


vector<double> moins(vector<double> u, vector<double> v){

    int n = v.size();
    int m = u.size();
    vector<double> res;
    res.resize(n);

    if (m==n){
        for(int i = 0; i<n;i ++){
            res[i] = u[i] - v[i];
        }
        return res;
    }
    else{
        cout << "les vecteurs u et v n'ont pas la mêmes dimensions pour le -"<< endl;
        return res;
    }

}






#define OPERATIONS_CPP
#endif
