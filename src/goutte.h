#ifndef GOUTTE_H

#include <iostream>
#include <vector>
#include <string>
#include "timescheme.h"


class goutte{

    private : 
        
        //Composante goutte
        int _id; 
        std::vector<double> _sol; 
        double _rho;
        double _d;  
        double _dt; 
        double _t0; 
        double _t; 
        double _tmax; 

        std::vector<double> _pos_imp; 
        std::vector<double> _pos_init; 

        //Composante air
        double _mu_a; 
        double _rho_a; 
        std::vector<double> _ua; 

        char _choix_interpolation; 
        std::string fichier_mesh; 
        std::string fichier_obstacle; 


    public : 

        goutte(int id); //Constructeur
        void initialisation(double d_min); 
        void advance(std::vector<std::vector<int>> tab_obstacle, std::vector<std::vector<double>> tab_noeud, std::vector<std::vector<int>> tab_triangle, std::vector<double> tab_rho_a, std::vector<std::vector<double>> tab_vitesse, std::vector<std::vector<int>> tab_voisin, double area_min, double vitesse_max); 



        double get_rho() const {return _rho;};
        double get_d() const {return _d;};
        double get_dt() const {return _dt;};
        double get_t0() const {return _t0;};
        double get_tmax() const {return _tmax;};



        double get_rho_a() const {return _rho_a;}; 
        double get_mu_a() const {return _mu_a;}; 
        std::vector<double> get_ua() const {return _ua;}; 
        std::vector<double> get_sol() const {return _sol;}; 
        std::vector<double> get_pos_imp() const {return _pos_imp;};
        std::vector<double> get_pos_init() const {return _pos_init;};

        std::string get_fichier() const {return fichier_mesh;}; 
        std::string get_fichier_obstacle() const {return fichier_obstacle;}; 

}; 


#define METHODE_H
#endif