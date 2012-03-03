/*
=================================================================
header.h
=================================================================
This is the head of the program, here the global variables 
are declared. The program was developed with three different 
integrators: RK4, RK-Felberg and RK-Cash-Karp. This methods 
are implemented in rk4.c and rk45.c files.

    Copyright (C) 2012  Edison Montoya, eamonto@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Up to date: 3 Mar 2012					
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define C        1.0                     //Speed of Light 
#define G        1.0 //6.6732D-11        //Newton's Constant
#define hbar     1.0 //1.05459D-34       //Planck's Constant
#define gamma    0.23753295796592        //Barbero-Immirzi Parameter

#define epsilon -1.0                     //Orientation Parameter
#define L1       1.0                     //Fiducial Lenght 1
#define L2       1.0                     //Fiducial Lenght 2
#define L3       1.0                     //Fiducial Lenght 3


///////////////// GLOBAL VARIABLES /////////////////

////// DYNAMICAL VARIABLES
double c1,c2,c3,p1,p2,p3;
double phi,P_phi;

////// IN-OUT VARIABLES
char *initfile,*outputfile;
int std_out, time_output;

////// TIME PARAMETERS
double initial_time, final_time, run_time, dt;
int  preferred_time;

////// PARAMETERS FOR THE DYNAMICAL EQUATIONS
int Integrator;
int dyn_eq;
double alpha;
double tol_err;
int itmax;
int pot_switch;
int pot_select;
double m_phi,l_phi;
double V0_pot,sigma1_pot,sigma2_pot;

////// OBSERVABLES  
double volume;
double H1, H2, H3;
double a1, a2, a3;
double a_prom;
double sigma2;
double omega;
double density;
double constraint;
double expansion;
double shear;
double Ricci;
double k1,k2,k3;
double curvature_param;

//AUXILIARY CONSTANTS
double Lp;
double V0;
double Delta;
double lambda;
double density_crit;
double onegamma2;
double invgamma; 
double initial_volume;

//AUXILIARY DYNAMICAL VARIABLES
double mu1, mu2, mu3;
double mu1c1, mu2c2, mu3c3;

////// AUXILIARIES FOR INTEGRATION 
double c1_p, c2_p, c3_p;
double p1_p, p2_p, p3_p;

double c1_a, c2_a, c3_a;
double p1_a, p2_a, p3_a;

double sc1, sc2, sc3;
double sp1, sp2, sp3;

double phi_p, phi_a, sphi;
double P_phi_p, P_phi_a, sP_phi;

//Adaptive RK
double c21;
double c31,c32;
double c41,c42,c43;
double c51,c52,c53,c54;
double c61,c62,c63,c64,c65;
double aa1,aa2,aa3,aa4,aa5,aa6;
double b1,b2,b3,b4,b5,b6;

//RKF Auxiliaries 
double c1_x4, c2_x4, c3_x4;
double p1_x4, p2_x4, p3_x4;

double c1_x5, c2_x5, c3_x5;
double p1_x5, p2_x5, p3_x5;

double phi_x4, phi_x5;
double P_phi_x4, P_phi_x5;

//Error
double max_err;

//RKF K's
double kc1[7],kc2[7],kc3[7];
double kp1[7],kp2[7],kp3[7];
double kphi[7];
double kP_phi[7];


///////////////LIBRARIES////////////////////////
#include <io_lib.h>
#include <initialize.h>
#include <observables.h>
#include <rk4.h>
#include <rk45.h>
#include <sources.h>

