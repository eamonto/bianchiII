/*
=========================================================================
rk45.c
=========================================================================
The Adaptive Runge-Kutta methods are implemented here. The method
is a RK45, which means that there are two Runge-Kutta, one of order 4 and 
other one of order 5. These kind of methods are the Runge-Kutta-Fehlberg 
and the Runge-Kutta-Cash-Karp.

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

///// STORE THE CONDITIONS BEFORE THE INTEGRATION
int store_levels_rk45(void)
{
  c1_p = c1;
  c2_p = c2;
  c3_p = c3;
  p1_p = p1;
  p2_p = p2;
  p3_p = p3;

  if(pot_switch==1)
    {
      phi_p   = phi;
      P_phi_p = P_phi;
    } 
  
  return 0;
}


///// FINAL VALUES AFTER INTEGRATION
int step_rk45(int i)
{
  if(i==0) //Repetion of iteration
    {
      c1 = c1_p;
      c2 = c2_p;
      c3 = c3_p;
      
      p1 = p1_p;
      p2 = p2_p;
      p3 = p3_p;

      if(pot_switch==1)
	{
	  phi   = phi_p;
	  P_phi = P_phi_p;
	} 
    }  
  else if(i==1) //Advanced in time
    {
      c1 = c1_x5;
      c2 = c2_x5;
      c3 = c3_x5;
      
      p1 = p1_x5;
      p2 = p2_x5;
      p3 = p3_x5;

      if(pot_switch==1)
	{
	  phi   = phi_x5;
	  P_phi = P_phi_x5;
	} 
    }

  return 0;
}


///// EVOLUTION WITH THE ADAPTIVE METHOD
int evolution_rk45(int i)
{
  kc1[i] = dt*sc1;
  kc2[i] = dt*sc2;
  kc3[i] = dt*sc3;
  
  kp1[i] = dt*sp1;
  kp2[i] = dt*sp2;
  kp3[i] = dt*sp3;
  
  if(pot_switch==1)
    {
      kphi[i]   = dt*sphi;
      kP_phi[i] = dt*sP_phi;
    }

  if(i==1)
    {
      c1 = c1_p + c21*kc1[i];
      c2 = c2_p + c21*kc2[i];
      c3 = c3_p + c21*kc3[i];

      p1 = p1_p + c21*kp1[i];
      p2 = p2_p + c21*kp2[i];
      p3 = p3_p + c21*kp3[i];

      if(pot_switch==1)
	{
	  phi   = phi_p   + c21*kphi[i];
	  P_phi = P_phi_p + c21*kP_phi[i];
	} 
    }
  else if(i==2)
    {
      c1 = c1_p + c31*kc1[1] +c32*kc1[2];
      c2 = c2_p + c31*kc2[1] +c32*kc2[2];
      c3 = c3_p + c31*kc3[1] +c32*kc3[2];
      
      p1 = p1_p + c31*kp1[1] +c32*kp1[2];
      p2 = p2_p + c31*kp2[1] +c32*kp2[2];
      p3 = p3_p + c31*kp3[1] +c32*kp3[2];

      if(pot_switch==1)
	{
	  phi   = phi_p   +c31*kphi[1]   +c32*kphi[2];
	  P_phi = P_phi_p +c31*kP_phi[1] +c32*kP_phi[2];
	} 
    }
  else if(i==3)
    {
      c1 = c1_p + c41*kc1[1] +c42*kc1[2] +c43*kc1[3];
      c2 = c2_p + c41*kc2[1] +c42*kc2[2] +c43*kc2[3];
      c3 = c3_p + c41*kc3[1] +c42*kc3[2] +c43*kc3[3];
      
      p1 = p1_p + c41*kp1[1] +c42*kp1[2] +c43*kp1[3];
      p2 = p2_p + c41*kp2[1] +c42*kp2[2] +c43*kp2[3];
      p3 = p3_p + c41*kp3[1] +c42*kp3[2] +c43*kp3[3];

      if(pot_switch==1)
	{
	  phi   = phi_p   +c41*kphi[1]   +c42*kphi[2]   +c43*kphi[3];
	  P_phi = P_phi_p +c41*kP_phi[1] +c42*kP_phi[2] +c43*kP_phi[3];
	} 
    }
  else if(i==4)
    {
      c1 = c1_p + c51*kc1[1] +c52*kc1[2] +c53*kc1[3] +c54*kc1[4];
      c2 = c2_p + c51*kc2[1] +c52*kc2[2] +c53*kc2[3] +c54*kc2[4];
      c3 = c3_p + c51*kc3[1] +c52*kc3[2] +c53*kc3[3] +c54*kc3[4];
      
      p1 = p1_p + c51*kp1[1] +c52*kp1[2] +c53*kp1[3] +c54*kp1[4];
      p2 = p2_p + c51*kp2[1] +c52*kp2[2] +c53*kp2[3] +c54*kp2[4];
      p3 = p3_p + c51*kp3[1] +c52*kp3[2] +c53*kp3[3] +c54*kp3[4];

      if(pot_switch==1)
	{
	  phi  = phi_p  +c51*kphi[1]  +c52*kphi[2]  +c53*kphi[3]  +c54*kphi[4];
	  P_phi= P_phi_p+c51*kP_phi[1]+c52*kP_phi[2]+c53*kP_phi[3]+c54*kP_phi[4];
	} 
    }
  else if(i==5)
    {
      c1 = c1_p +c61*kc1[1] +c62*kc1[2] +c63*kc1[3] +c64*kc1[4] +c65*kc1[5];
      c2 = c2_p +c61*kc2[1] +c62*kc2[2] +c63*kc2[3] +c64*kc2[4] +c65*kc2[5];
      c3 = c3_p +c61*kc3[1] +c62*kc3[2] +c63*kc3[3] +c64*kc3[4] +c65*kc3[5];
      
      p1 = p1_p +c61*kp1[1] +c62*kp1[2] +c63*kp1[3] +c64*kp1[4] +c65*kp1[5];
      p2 = p2_p +c61*kp2[1] +c62*kp2[2] +c63*kp2[3] +c64*kp2[4] +c65*kp2[5];
      p3 = p3_p +c61*kp3[1] +c62*kp3[2] +c63*kp3[3] +c64*kp3[4] +c65*kp3[5];

      if(pot_switch==1)
	{
	  phi  = phi_p   +c61*kphi[1]   +c62*kphi[2]   +c63*kphi[3]   +c64*kphi[4]   +c65*kphi[5];
	  P_phi= P_phi_p +c61*kP_phi[1] +c62*kP_phi[2] +c63*kP_phi[3] +c64*kP_phi[4] +c65*kP_phi[5];
	} 
    }

  return 0;
}


///// RELATIVE ERROR BETWEEN THE RK4 AND RK5 METHODS
int error_rk45(void)
{
  double aux_err[8];
  int i;

  for(i=0;i<8;i++) aux_err[i]=0.0;

  c1_x4 = c1_p +aa1*kc1[1] +aa3*kc1[3] +aa4*kc1[4] +aa5*kc1[5] +aa6*kc1[6];
  c2_x4 = c2_p +aa1*kc2[1] +aa3*kc2[3] +aa4*kc2[4] +aa5*kc2[5] +aa6*kc2[6];
  c3_x4 = c3_p +aa1*kc3[1] +aa3*kc3[3] +aa4*kc3[4] +aa5*kc3[5] +aa6*kc3[6];
  
  p1_x4 = p1_p +aa1*kp1[1] +aa3*kp1[3] +aa4*kp1[4] +aa5*kp1[5] +aa6*kp1[6];
  p2_x4 = p2_p +aa1*kp2[1] +aa3*kp2[3] +aa4*kp2[4] +aa5*kp2[5] +aa6*kp2[6];
  p3_x4 = p3_p +aa1*kp3[1] +aa3*kp3[3] +aa4*kp3[4] +aa5*kp3[5] +aa6*kp3[6];
  
  
  c1_x5 = c1_p +b1*kc1[1] +b3*kc1[3] +b4*kc1[4] +b5*kc1[5] +b6*kc1[6];
  c2_x5 = c2_p +b1*kc2[1] +b3*kc2[3] +b4*kc2[4] +b5*kc2[5] +b6*kc2[6];
  c3_x5 = c3_p +b1*kc3[1] +b3*kc3[3] +b4*kc3[4] +b5*kc3[5] +b6*kc3[6];
  
  p1_x5 = p1_p +b1*kp1[1] +b3*kp1[3] +b4*kp1[4] +b5*kp1[5] +b6*kp1[6];
  p2_x5 = p2_p +b1*kp2[1] +b3*kp2[3] +b4*kp2[4] +b5*kp2[5] +b6*kp2[6];
  p3_x5 = p3_p +b1*kp3[1] +b3*kp3[3] +b4*kp3[4] +b5*kp3[5] +b6*kp3[6];

  aux_err[0] = fabs(c1_x4-c1_x5);
  aux_err[1] = fabs(c2_x4-c2_x5);
  aux_err[2] = fabs(c3_x4-c3_x5);
  aux_err[3] = fabs(p1_x4-p1_x5);
  aux_err[4] = fabs(p2_x4-p2_x5);
  aux_err[5] = fabs(p3_x4-p3_x5);
  
  if(pot_switch==1)
    {
      phi_x4   = phi_p   +aa1*kphi[1]   +aa3*kphi[3]   +aa4*kphi[4]   +aa5*kphi[5]   +aa6*kphi[6];
      P_phi_x4 = P_phi_p +aa1*kP_phi[1] +aa3*kP_phi[3] +aa4*kP_phi[4] +aa5*kP_phi[5] +aa6*kP_phi[6];
      
      phi_x5   = phi_p   +b1*kphi[1]   +b3*kphi[3]   +b4*kphi[4]   +b5*kphi[5]   +b6*kphi[6];
      P_phi_x5 = P_phi_p +b1*kP_phi[1] +b3*kP_phi[3] +b4*kP_phi[4] +b5*kP_phi[5] +b6*kP_phi[6];
      
      aux_err[6] = fabs(phi_x4-phi_x5);
      aux_err[7] = fabs(P_phi_x4-P_phi_x5);
    }

  max_err = 0.0;

  for(i=0;i<8;i++)
    {
      if(aux_err[i] > max_err)
	max_err = aux_err[i];
    }

  return 0;
}

