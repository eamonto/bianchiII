/*
========================================================================
main.c
========================================================================
Principal program. This program solves the Bianchi I and II equations,
classical and effective ones, without potential or with inflationary 
and cyclic potentials.

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


int main (int argc,char **argv)
{
  double Ttotal,s;
  int k,l,m; 
  char newfile[200],*aux;  
  FILE *pf;

  k=0,
  l=0;
  m=0;
  s=0.0;

  if(argc != 3) usage();        //VERIFICATION OF INPUT FILES

  initfile=argv[1];             //ASIGNATION OF FILE'S 
  outputfile=argv[2];           //NAMES
  
  read_param();                 //READ INITIAL CONDITIONS

  create_remove_dir();          //REMOVE OLD outputfile AND CREATE A NEW ONE

  initialize_all();             //INITIALIZATION OF ALL QUANTITIES

  check_initial_data();         //CHECK IF THE INITIAL CONDITIONS ARE CORRECT

  compute_obs();                //COMPUTE OBSERVABLES

  write_output();               //WRITE THE OBSERVABLES TO A FILE

  aux="/rk45.dat";              //FILE FOR INFORMATION
  strcpy(newfile,outputfile);   //OF THE ADAPTIVE
  strcat(newfile,aux);          //INTEGRATORS

  if(std_out==1){
    printf("\n \t Start Simulation \n");
    printf("\t Time -> %4.6lf\n",run_time);
  }

  Ttotal = final_time-initial_time;

  do
    {      
      if(Integrator==0)             //RK4 INTEGRATOR
	{
	  store_levels_rk4();

	  rhs();
	  evolution_rk4(1);
	  rhs();
	  evolution_rk4(2);
	  rhs();
	  evolution_rk4(3);
	  rhs();
	  evolution_rk4(4);
	  
	  run_time = run_time + dt;	  
	}
      
      else if(Integrator > 0)      //ADAPTIVE INTEGRATOR RK45
	{
	  store_levels_rk45();
	  
	  for(k=1;k<=6;k++){
	    rhs();
	    evolution_rk45(k);
	  }
	  
	  error_rk45();            // RELATIVE ERROR
	  
	  if(max_err==0.0) 
	    {
	      if(std_out==1)
		printf("\t The error is zero at iteration %d\n",l);
	      
	      step_rk45(1);        //EVOLVE THE SYSTEM
	      run_time = run_time + dt;
	      dt = 2.0*dt;
	    }
	  
	  else
	    {
	      if(max_err > tol_err)
		{
		  step_rk45(0);    //START AGAIN
		  s = 0.84*pow(fabs(tol_err/max_err),0.25);
		  dt= dt*s;        //PREDICTION OF THE NEW TIME STEP
		  m++;
		}
	      else
		{
		  step_rk45(1);     //EVOLVE THE SYSTEM
		  run_time = run_time + dt;
		  s = 0.84*pow(fabs(tol_err/max_err),0.2);
		  dt= dt*s;        //PREDICTION OF THE NEW TIME STEP
		  m = 0;
		}
	    }
	}
      
      l++;     

      if(m>itmax)
	{
	  printf("\t Too much iterations at time %4.6lf\n",run_time);
	  exit(0);
	}
      
      if(fabs(dt)<1.0e-10)
	{
	  printf("\t Too small dt=%e\n",dt);
	  exit(0);
	}
      
      if(l%time_output==0)
	{
	  compute_obs();          //COMPUTE OBSERVABLES
	  write_output();         //WRITE THE OBSERVABLES TO A FILE
	  if(std_out==1)
	    printf("\t Time -> %4.6lf\n",run_time);
	}
      
      //INFORMATION OF THE ADAPTIVE INTEGRATOR 
      if(Integrator > 0)
	{
	  pf=fopen(newfile,"a");
	  fprintf(pf,"%.15e\t%.15e\t%.15e\t%.15e\t%d\t%d\n",
			 run_time,fabs(dt),s,max_err,m,l);
	  fclose(pf);
	}
      
    }while(fabs(Ttotal)>fabs(initial_time-run_time));

  if(std_out==1)  printf("\n \t Finish Simulation \n");
  
  return 0;
}
