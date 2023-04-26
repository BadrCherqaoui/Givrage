#ifndef TIMESCHEME_CPP

#include "goutte.h"
#include "timescheme.h"
#include "operations.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

//=========================== Calcul de Re et Cd défini dans la thèse =============================


double calcul_Re(vector<double> sol,double d, double rho, double rho_a, double mu_a, vector<double> ua){
    
    double norme;
    int dim(ua.size());
    vector<double> U(dim, 0);

    for(int i=dim; i<sol.size(); i++) { U[i-dim]=sol[i]; }

    norme = norm(moins(ua, U));
    double Re = (rho_a*norme*d)/mu_a; 

    return Re; 
}

double calcul_Cd(double Re) {
    double Cd;

    if (Re < 1e-11){
        Cd = 0; 
    }
    else if (Re <= 0.5){
        Cd = 24./Re; 
    }
    else if (Re <= 10){
        Cd = 27.17/pow(Re, 0.821);
    }
    else if (Re <= 100){
        Cd = 15.268/pow(Re, 0.571);
    }
    else if (Re <= 1000){
        Cd = 6.877/pow(Re, 0.398);
    }
    else {
        Cd = 0.44; 
    }

    return Cd;

}



// ==========================  Création de f de l'équation différentiel y' = f(t,y) =============================

vector<double> f(double t,vector<double> sol, double d, double rho, double rho_a, double mu_a, vector<double> ua, double Re, double Cd){

    int dim(ua.size()); 
    vector<double> res(2*dim, 0); 
    //Calcul de f 

    for (int i=0; i<dim; i++) {
        res[i] = sol[i + dim];
        res[i + dim] = 3.*mu_a*Cd*Re/(rho*d*d*4.)*(ua[i]-sol[i + dim]);
    }

    return res; 
}



// ========================= RK2 et RK4 pour résoudre le problème lagrangien =========================

void RK2(double t,vector<double> &sol, double dt, double d, double rho, double rho_a, double mu_a, vector<double> ua, double & Re)
{  
    int dim_sol(sol.size());

    double Cd; 
    Re = calcul_Re(sol, d,rho, rho_a, mu_a,ua); 
    Cd = calcul_Cd(Re); 
    
    vector<double> k(dim_sol,0), kinter(dim_sol,0), kinter2(dim_sol,0);
    kinter = prod_const_vector(dt/2., f(t,sol,d,rho, rho_a, mu_a, ua, Re, Cd));
    kinter2 = addition(sol, kinter);
    k = f(t+dt/2., kinter2,d,rho, rho_a, mu_a, ua, Re, Cd);

    sol = addition(sol, prod_const_vector(dt, k));

}

void RK4(double t,vector<double> &sol, double dt, double d, double rho, double rho_a, double mu_a, vector<double> ua, double & Re){   

    int dim_sol(sol.size());
    vector<double> k1(dim_sol,0),k2(dim_sol,0),k3(dim_sol,0),k4(dim_sol,0),ktot(dim_sol,0);

    double Cd; 
    Re = calcul_Re(sol, d, rho, rho_a, mu_a, ua); 
    Cd = calcul_Cd(Re); 


    k1=f(t,sol,d,rho, rho_a, mu_a, ua,Re, Cd);
    ktot=k1;
    k2=f(t+dt/2., addition(sol,prod_const_vector(dt/2.,k1)),d,rho, rho_a, mu_a, ua, Re, Cd);
    ktot=addition(prod_const_vector(2.,k2),ktot);
    k3=f(t+dt/2., addition(sol,prod_const_vector(dt/2.,k2)),d,rho, rho_a, mu_a, ua, Re, Cd);
    ktot=addition(prod_const_vector(2.,k3),ktot);
    k4=f(t+dt, addition(sol,prod_const_vector(dt,k3)),d,rho, rho_a, mu_a, ua, Re, Cd);
    ktot=addition(k4,ktot);

    sol=addition(sol,prod_const_vector(dt/6.,ktot));


}    








#define TIMESCHEME_CPP
#endif
