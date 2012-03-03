/*
=========================================================================
initialize.c
=========================================================================
Initialization of all quantities, including the parameters for
the adaptive methods (Fehlberg and Cash-Karp).

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

#include <header.h>


////// INITIALIZATION OF ALL QUANTITIES
int initialize_all(void)
{
  int i;
  char newfile[200],*aux;  
  FILE *pf;

  //Observables
  volume = 1.0;
  H1 = 0.0;
  H2 = 0.0;
  H3 = 0.0;
  a1 = 0.0;
  a2 = 0.0;
  a3 = 0.0;
  a_prom = 0.0;
  sigma2 = 0.0;
  omega = 0.0;
  density = 0.0;
  constraint = 0.0;
  expansion = 0.0;
  shear = 0.0;
  Ricci = 0.0;
  k1 = 0.0;
  k2 = 0.0;
  k3 = 0.0;
  curvature_param = 0.0;

  //Time
  run_time = initial_time;

  //Auxiliary constants
  Lp = sqrt(hbar*G/(C*C*C)); 
  V0 = L1*L2*L3;                
  Delta = 4.0*M_PI*sqrt(3.0)*gamma*Lp*Lp; //5.170045537718 
  lambda = sqrt(Delta);                   //2.273773413891       
  density_crit = sqrt(3.0)/(32.0*M_PI*M_PI*gamma*gamma*gamma*G*G*hbar);
  onegamma2 = 1.0+gamma*gamma;
  invgamma = 1.0/gamma;  
  initial_volume = sqrt(p1*p2*p3);

  //Step
  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;

  //Dynamical variables
  c1 = mu1c1/mu1;
  c2 = mu2c2/mu2;
  c3 = mu3c3/mu3;

  if(std_out==1){
    printf("\n c1 -> %lf",c1);
    printf("\n c2 -> %lf",c2);
    printf("\n c3 -> %lf",c3);
    printf("\n ");    
  }
  
  //Initial Field Momentum
  P_phi = field_momentum();  
  if(std_out==1) printf("\n Field Momentum -> %lf \n",P_phi);

  //Auxiliary for integration
  c1_p = 0.0;
  c2_p = 0.0;  
  c3_p = 0.0;
  p1_p = 0.0;
  p2_p = 0.0;
  p3_p = 0.0;
  phi_p = 0.0;
  P_phi_p = 0.0;
  
  c1_a = 0.0;
  c2_a = 0.0; 
  c3_a = 0.0;
  p1_a = 0.0; 
  p2_a = 0.0; 
  p3_a = 0.0;
  phi_a = 0.0;
  P_phi_a = 0.0;
  
  sc1 = 0.0;
  sc2 = 0.0;
  sc3 = 0.0;
  sp1 = 0.0;
  sp2 = 0.0;
  sp3 = 0.0;
  sphi = 0.0;
  sP_phi = 0.0;
  
  //Adaptive RK
  if (Integrator==1) {  //Fehlberg
    
    c21 = 1.0/4.0;
    c31 = 3.0/32.0;
    c32 = 9.0/32.0;
    c41 = 1932.0/2197.0;
    c42 = -7200.0/2197.0;
    c43 = 7296.0/2197.0;
    c51 = 439.0/216.0;
    c52 = -8.0;
    c53 = 3680.0/513.0;
    c54 = -845.0/4104.0;
    c61 = -8.0/27.0;
    c62 = 2.0;
    c63 = -3544.0/2565.0;
    c64 = 1859.0/4104.0;
    c65 = -11.0/40.0;
    
    aa1 = 25.0/216.0;
    aa2 = 0.0;
    aa3 = 1408.0/2565.0;
    aa4 = 2197.0/4104.0;
    aa5 = -1.0/5.0;
    aa6 = 0.0;
    
    b1 = 16.0/135.0;
    b2 = 0.0;
    b3 = 6656.0/12825.0;
    b4 = 28561.0/56430.0;
    b5 = -9.0/50.0;
    b6 = 2.0/55.0;
    
  }else if (Integrator==2) { //CASH-KARP
    
    c21 = 1.0/5.0;
    c31 = 3.0/40.0;
    c32 = 9.0/40.0;
    c41 = 3.0/10.0;
    c42 = -9.0/10.0;
    c43 = 6.0/5.0;
    c51 = -11.0/54.0;
    c52 = 5.0/2.0;
    c53 = -70.0/27.0;
    c54 = 35.0/27.0;
    c61 = 1631.0/55296.0;
    c62 = 175.0/512.0;
    c63 = 575.0/13824.0;
    c64 = 44275.0/110592.0;
    c65 = 253.0/4096.0;
    
    aa1 = 37.0/378.0;
    aa2 = 0.0;
    aa3 = 250.0/621.0;
    aa4 = 125.0/594.0;
    aa5 = 0.0;
    aa6 = 512.0/1771.0;
    
    b1 = 2825.0/27648.0;
    b2 = 0.0;
    b3 = 18575.0/48384.0;
    b4 = 13525.0/55296.0;
    b5 = 277.0/14336.0;
    b6 = 1.0/4.0;
    
  }
  
  //RKF Auxiliaries 
  c1_x4 = 0.0;
  c2_x4 = 0.0; 
  c3_x4 = 0.0;
  p1_x4 = 0.0;
  p2_x4 = 0.0; 
  p3_x4 = 0.0;
  phi_x4 = 0.0;
  P_phi_x4 = 0.0;
  
  c1_x5 = 0.0; 
  c2_x5 = 0.0; 
  c3_x5 = 0.0;
  p1_x5 = 0.0; 
  p2_x5 = 0.0; 
  p3_x5 = 0.0;
  phi_x5 = 0.0;
  P_phi_x5 = 0.0;

  //RKF K's
  for(i=0;i<7;i++)
    {
      kc1[i] = 0.0;
      kc2[i] = 0.0;
      kc3[i] = 0.0;
      kp1[i] = 0.0;
      kp2[i] = 0.0;
      kp3[i] = 0.0;
      kphi[i] = 0.0;
      kP_phi[i] = 0.0;
    }

  return 0;
}
