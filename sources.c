/*
=========================================================================
sources.c
=========================================================================
Here are implemented the right hand side of the equations 
that we want to integrate.

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


///// RIGHT HAND SIDE /////
int rhs(void)
{
  /////RHS OF THE CLASSICAL EQUATIONS////
  if(dyn_eq==0)
    {
      //Internal Time
      
      sc1 = -invgamma*( p2*c1*c2 + p3*c1*c3 
			+alpha*alpha*0.5*onegamma2*(p2*p3/p1)*(p2*p3/p1)/p1 )
	    +8.0*M_PI*G*gamma*p2*p3*potential();
      
      sc2 = -invgamma*( p1*c2*c1 + p3*c2*c3 + alpha*epsilon*p3*c1  
			-alpha*alpha*0.5*onegamma2*(p3/p1)*(p3/p1)*p2 )
	    +8.0*M_PI*G*gamma*p1*p3*potential();
      
      sc3 = -invgamma*( p1*c3*c1 + p2*c3*c2 + alpha*epsilon*p2*c1 
			-alpha*alpha*0.5*onegamma2*(p2/p1)*(p2/p1)*p3 )
	    +8.0*M_PI*G*gamma*p1*p2*potential();
      
          
      sp1 = invgamma*(p1*p2*c2 + p1*p3*c3 + alpha*epsilon*p2*p3);
      
      sp2 = invgamma*(p2*p1*c1 + p2*p3*c3);
      
      sp3 = invgamma*(p3*p1*c1 + p3*p2*c2);

      if(pot_switch==1)
	{
	  sphi   = P_phi;
	  sP_phi =-p1*p2*p3*dev_pot();
	}


      if(preferred_time==1)//Cosmic Time
	{
	  volume = sqrt(p1*p2*p3);
      
	  sc1 = sc1/volume;
	  sc2 = sc2/volume;
	  sc3 = sc3/volume;

	  sp1 = sp1/volume;
	  sp2 = sp2/volume;
	  sp3 = sp3/volume;

	  if(pot_switch==1)
	    {
	      sphi   = sphi/volume;
	      sP_phi = sP_phi/volume;
	    }
	}
    }  

  /////RHS OF THE EFFECTIVE EQUATIONS////       
  else if(dyn_eq==1)
    {
      mu1 = sqrt(p1/(p2*p3))*lambda;
      mu2 = sqrt(p2/(p1*p3))*lambda;
      mu3 = sqrt(p3/(p1*p2))*lambda;
    
      //Internal Time
      
      sc1 = -invgamma*( p2*p3/(Delta)*
		        ( sin(mu1*c1)*sin(mu2*c2)    
		        + sin(mu1*c1)*sin(mu3*c3) 
		        + sin(mu2*c2)*sin(mu3*c3) 
		        + 0.5*mu1*c1*cos(mu1*c1)*(sin(mu2*c2)+sin(mu3*c3)) 
		        - 0.5*mu2*c2*cos(mu2*c2)*(sin(mu1*c1)+sin(mu3*c3)) 
		        - 0.5*mu3*c3*cos(mu3*c3)*(sin(mu1*c1)+sin(mu2*c2)) 
		        ) 
		        + 0.5*alpha*alpha*onegamma2*(p2*p3/p1)*(p2*p3/p1)/p1 
			+ 0.5*alpha/(lambda)*pow(p2*p3/p1,1.5) 
			*(mu1*c1*cos(mu1*c1)-sin(mu1*c1)) 
		      )
	    +8.0*M_PI*G*gamma*p2*p3*potential();
      
      sc2 = -invgamma*( p1*p3/(Delta)*
			( sin(mu1*c1)*sin(mu2*c2)  
		        + sin(mu1*c1)*sin(mu3*c3) 
		        + sin(mu2*c2)*sin(mu3*c3) 
		        - 0.5*mu1*c1*cos(mu1*c1)*(sin(mu2*c2)+sin(mu3*c3)) 
		        + 0.5*mu2*c2*cos(mu2*c2)*(sin(mu1*c1)+sin(mu3*c3)) 
			- 0.5*mu3*c3*cos(mu3*c3)*(sin(mu1*c1)+sin(mu2*c2)) 
			)
			- 0.5*alpha*alpha*onegamma2*p2*p3*p3/(p1*p1) 
			+ 0.5*alpha*p3/mu1 
			*(3.0*sin(mu1*c1)-mu1*c1*cos(mu1*c1)) 
		      )
	    +8.0*M_PI*G*gamma*p1*p3*potential();
      
      sc3 = -invgamma*( p1*p2/(Delta)*
			( sin(mu1*c1)*sin(mu2*c2)  
			+ sin(mu1*c1)*sin(mu3*c3) 
			+ sin(mu2*c2)*sin(mu3*c3) 
			- 0.5*mu1*c1*cos(mu1*c1)*(sin(mu2*c2)+sin(mu3*c3)) 
			- 0.5*mu2*c2*cos(mu2*c2)*(sin(mu1*c1)+sin(mu3*c3)) 
			+ 0.5*mu3*c3*cos(mu3*c3)*(sin(mu1*c1)+sin(mu2*c2)) 
			)
			- 0.5*alpha*alpha*onegamma2*p2*p2*p3/(p1*p1) 
			+ 0.5*alpha*p2/mu1 
			*(3.0*sin(mu1*c1)-mu1*c1*cos(mu1*c1)) 
		      )
	    +8.0*M_PI*G*gamma*p1*p2*potential();
      
      
      sp1 = invgamma*( p1*p1/mu1*(sin(mu2*c2) + sin(mu3*c3))+alpha*p2*p3 )*cos(mu1*c1);
      
      sp2 = invgamma*p2*p2*( sin(mu1*c1) + sin(mu3*c3) )*cos(mu2*c2)/mu2;
      
      sp3 = invgamma*p3*p3*( sin(mu1*c1) + sin(mu2*c2) )*cos(mu3*c3)/mu3;
      
      if(pot_switch==1)
	{
	  sphi   = P_phi;
	  sP_phi = -p1*p2*p3*dev_pot();
	}


      if(preferred_time==1)//Cosmic Time
	{
	  volume = sqrt(p1*p2*p3);
	  
	  sc1 = sc1/volume;
	  sc2 = sc2/volume;
	  sc3 = sc3/volume;
	  
	  sp1 = sp1/volume;
	  sp2 = sp2/volume;
	  sp3 = sp3/volume;

	  if(pot_switch==1)
	    {
	      sphi   = sphi/volume;
	      sP_phi = sP_phi/volume;
	    }
	}
    }
  
  return 0;
}
