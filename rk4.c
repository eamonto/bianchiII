/*
=========================================================================
rk4.c
=========================================================================
Implementation of the Runge-Kutta 4 method.

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
int store_levels_rk4(void)
{
  c1_p = c1;
  c2_p = c2;  
  c3_p = c3;
  p1_p = p1;
  p2_p = p2;
  p3_p = p3;
  
  c1_a = c1;
  c2_a = c2;  
  c3_a = c3;
  p1_a = p1;
  p2_a = p2;
  p3_a = p3;

  if(pot_switch==1)
    {
      phi_p   = phi;
      P_phi_p = P_phi;

      phi_a   = phi;
      P_phi_a = P_phi;
    } 
  
  return 0;
}


///// EVOLUTION WITH THE RK4 METHOD
int evolution_rk4(int k)
{
  double dtw,weight;

  if (k==1)
    {
      dtw = dt*0.5;
      weight = dt/6.0;

      c1 = c1_p + dtw*sc1;
      c2 = c2_p + dtw*sc2;
      c3 = c3_p + dtw*sc3;
      
      p1 = p1_p + dtw*sp1;
      p2 = p2_p + dtw*sp2;
      p3 = p3_p + dtw*sp3;
      
      c1_a = c1_a + weight*sc1;
      c2_a = c2_a + weight*sc2;
      c3_a = c3_a + weight*sc3;
      
      p1_a = p1_a + weight*sp1;
      p2_a = p2_a + weight*sp2;
      p3_a = p3_a + weight*sp3;

      if(pot_switch==1)
	{
	  phi   = phi_p   + dtw*sphi;
	  P_phi = P_phi_p + dtw*sP_phi;

	  phi_a   = phi_a   + weight*sphi;
	  P_phi_a = P_phi_a + weight*sP_phi;
  	}
    }
  else if (k==2)
    {
      dtw = dt*0.5;
      weight = dt/3.0;

      c1 = c1_p + dtw*sc1;
      c2 = c2_p + dtw*sc2;
      c3 = c3_p + dtw*sc3;
      
      p1 = p1_p + dtw*sp1;
      p2 = p2_p + dtw*sp2;
      p3 = p3_p + dtw*sp3;
      
      c1_a = c1_a + weight*sc1;
      c2_a = c2_a + weight*sc2;
      c3_a = c3_a + weight*sc3;
      
      p1_a = p1_a + weight*sp1;
      p2_a = p2_a + weight*sp2;
      p3_a = p3_a + weight*sp3;

      if(pot_switch==1)
	{
	  phi   = phi_p   + dtw*sphi;
	  P_phi = P_phi_p + dtw*sP_phi;

	  phi_a   = phi_a   + weight*sphi;
	  P_phi_a = P_phi_a + weight*sP_phi;	  
	}
    }
  else if (k==3)
    {
      dtw = dt;
      weight = dt/3.0;

      c1 = c1_p + dtw*sc1;
      c2 = c2_p + dtw*sc2;
      c3 = c3_p + dtw*sc3;
  
      p1 = p1_p + dtw*sp1;
      p2 = p2_p + dtw*sp2;
      p3 = p3_p + dtw*sp3;
      
      c1_a = c1_a + weight*sc1;
      c2_a = c2_a + weight*sc2;
      c3_a = c3_a + weight*sc3;
      
      p1_a = p1_a + weight*sp1;
      p2_a = p2_a + weight*sp2;
      p3_a = p3_a + weight*sp3;

      if(pot_switch==1)
	{
	  phi   = phi_p   + dtw*sphi;
	  P_phi = P_phi_p + dtw*sP_phi;

	  phi_a   = phi_a   + weight*sphi;
	  P_phi_a = P_phi_a + weight*sP_phi;	  	  
	}
    }
  else if (k==4) 
    {
      weight = dt/6.0;
      
      c1 = c1_a + weight*sc1;
      c2 = c2_a + weight*sc2;
      c3 = c3_a + weight*sc3;
      
      p1 = p1_a + weight*sp1;
      p2 = p2_a + weight*sp2;
      p3 = p3_a + weight*sp3;

      if(pot_switch==1)
	{
	  phi   = phi_a   + weight*sphi;
	  P_phi = P_phi_a + weight*sP_phi;	  	  	  
	}
    }

  return 0;
}

