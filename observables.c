/*
=========================================================================
observables.c
=========================================================================
Here are compute all the quantities that can give relevant 
information about the dynamics of the system.

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

////// SCALAR FIELD POTENTIAL
double potential(void)
{
  double pot = 0.0;

  if(pot_switch==1)
    {
      if(pot_select==0){ //Inflationary potential
	pot = 0.5*m_phi*m_phi*phi*phi + l_phi*phi*phi*phi*phi/24.0;
      }else{             //Cyclic potential
	pot = V0_pot*(1.0 - exp(-sigma1_pot*phi))*exp( -exp(-sigma2_pot*phi) );
      }
    }

  return pot;
}


////// DERIVATIVE OF THE SCALAR FIELD POTENTIAL
double dev_pot(void)
{
  double pot = 0.0;

  if(pot_switch==1)
    {
      if(pot_select==0){ //Inflationary potential
	pot = m_phi*m_phi*phi + l_phi*phi*phi*phi/6.0;
      }else{             //Cyclic potential
	pot = potential()*( sigma1_pot/(exp(sigma1_pot*phi) -1.0) 
			   +sigma2_pot*exp(-sigma2_pot*phi)
			  ); 
      }
    }

  return pot;
}


////// SCALAR FIELD MOMENTUM
double field_momentum(void)
{
  double P_phi_aux = 0.0;

  //SQRT ARGUMENT FOR THE MOMENTUM
  if (dyn_eq==0) {       //Classical    
    
    P_phi_aux = 1.0/(8.0*M_PI*G*gamma*gamma)
                  *(p1*p2*c1*c2 + p2*p3*c2*c3 + p1*p3*c1*c3 
		    + alpha*epsilon*p2*p3*c1 
		    - (1.0+gamma*gamma)*(alpha*p2*p3/(2.0*p1))
		    *(alpha*p2*p3/(2.0*p1)) 
		    ) -p1*p2*p3*potential();
        
  }else if (dyn_eq==1) { //Effective
        
    P_phi_aux = p1*p2*p3/(8.0*M_PI*G*gamma*gamma*Delta) 
			     *( sin(mu1*c1)*sin(mu2*c2) 
			       +sin(mu2*c2)*sin(mu3*c3) 
			       +sin(mu3*c3)*sin(mu1*c1)
			      ) 
               + 1.0/(8.0*M_PI*G*gamma*gamma) 
		  *( alpha*pow(p2*p3,1.5)/(sqrt(Delta*p1))*sin(mu1*c1) 
		    -(1.0+gamma*gamma)*(alpha*p2*p3/(2.0*p1))
		     *(alpha*p2*p3/(2.0*p1)) 
		   ) -p1*p2*p3*potential();
  }

  if (P_phi_aux < 0.0){
    printf("\n\t WARNING, the sqrt argument is negative = %lf \n\n",P_phi_aux);
    exit(1);
  }

  P_phi_aux = sqrt(2.0*P_phi_aux);

  return P_phi_aux;
}


////// COMPUTE OBSERVABLES
int compute_obs(void)
{
  double thetadot=0.0,funx=0.0;
  double P_phi_aux;
  int time_aux;     
  
  mu1 = sqrt(p1/(p2*p3))*lambda;
  mu2 = sqrt(p2/(p1*p3))*lambda;
  mu3 = sqrt(p3/(p1*p2))*lambda;
  
  a1 = sqrt(p2*p3/p1)/L1;
  a2 = sqrt(p1*p3/p2)/L2;
  a3 = sqrt(p1*p2/p3)/L3;

  a_prom = pow(a1*a2*a3,1.0/3.0);
  
  volume = sqrt(p1*p2*p3);
  
  time_aux = preferred_time;
  preferred_time = 1; //The observables are computed in the cosmic time (Lapse=1)
  rhs();  //Sources in the cosmic time
  preferred_time = time_aux;
  
  H1 = (sp2/p2 + sp3/p3 - sp1/p1)/2.0;
  H2 = (sp1/p1 + sp3/p3 - sp2/p2)/2.0;
  H3 = (sp1/p1 + sp2/p2 - sp3/p3)/2.0;
  
  expansion = H1+H2+H3;
  
  shear = ( (H1-H2)*(H1-H2) + (H2-H3)*(H2-H3) + (H3-H1)*(H3-H1) )/3.0;

  P_phi_aux = field_momentum();
  
  constraint = fabs(P_phi_aux-P_phi)/P_phi_aux;
  
  density = P_phi*P_phi/(2.0*p1*p2*p3) + potential(); 

  if(fabs(expansion) > 1.0e-15) //Avoids division by zero
    {
      omega = 24.0*M_PI*G*density/(expansion*expansion);
      
      sigma2 = 3.0*shear/(2.0*expansion*expansion);

      curvature_param = 3.0*p2*p3/(4.0*p1*p1*p1*expansion*expansion);

      // Kasner Exponents
      k1 = H1/fabs(expansion);
      k2 = H2/fabs(expansion);
      k3 = H3/fabs(expansion);
    }
  
  if (dyn_eq==0){ //CLASSICAL

    //Raychaudhuri Equation
    thetadot = -0.5*expansion*expansion -shear -16.0*M_PI*G*density; 

    if (density/density_crit > 10.0){
      printf("\n\t We are near to the Big Bang.");
      printf("\n\t The program stops.\n");
      exit(0);
    }

  }else if(dyn_eq==1){ //EFFECTIVE
    double mu1dot,mu2dot,mu3dot;
    double sin1dot,sin2dot,sin3dot;
    
    mu1dot = -mu1*H1;
    mu2dot = -mu2*H2;
    mu3dot = -mu3*H3;
    
    sin1dot = cos(mu1*c1)*(mu1dot*c1 + mu1*sc1);
    sin2dot = cos(mu2*c2)*(mu2dot*c2 + mu2*sc2);
    sin3dot = cos(mu3*c3)*(mu3dot*c3 + mu3*sc3);

    funx = alpha*sqrt(p2*p3/(p1*p1*p1));
        
    thetadot = 1.0/(2.0*gamma*lambda)*
                ( 
		  2.0*sin(mu1*c1) + cos(mu1*c1)*(sin2dot+sin3dot)
		 +2.0*sin(mu2*c2) + cos(mu2*c2)*(sin1dot+sin3dot)
		 +2.0*sin(mu3*c3) + cos(mu3*c3)*(sin1dot+sin2dot)
		 +lambda*funx*(1.0 + 0.5*cos(mu1*c1)*(sp2/p2+sp3/p3-3.0*sp1/p1) )   
		);
  }
   
  Ricci = 2.0*thetadot +(sp1/p1)*(sp1/p1) +(sp2/p2)*(sp2/p2) +(sp3/p3)*(sp3/p3)-0.5*funx*funx;

  //  Ricci = 2.0*thetadot +shear +4.0*expansion*expansion/3.0 -0.5*funx*funx;

  return 0;
}
