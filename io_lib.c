/*
=========================================================================
io_lib.c
=========================================================================

This contains:
-> Input-output routines.
-> Create and remove files routines.
-> Usage message routine.

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


//REMOVE OLD outputfile AND CREATE A NEW ONE
int create_remove_dir(void)
{
  char remove[200],copy[200],create[200];
  char *rm="rm -rf ",*cp="cp ",*mkdir="mkdir -p ";
  int out_return;

  strcpy(remove,rm);
  strcat(remove,outputfile);
  out_return = system(remove);

  strcpy(create,mkdir);
  strcat(create,outputfile);
  out_return = system(create);

  strcpy(copy,cp);
  strcat(copy,initfile);
  strcat(copy," ");
  strcat(copy,outputfile);
  out_return = system(copy);

  return out_return;
}


///////////VERIFICATION OF INPUT FILES/////////////
int usage(void)
{
  printf("\n");
  printf("USAGE: ./exec <initfile> <outputfile> \n");
  printf("\n");
  exit(1);
}


///////////CHECK IF THE INITIAL CONDITIONS ARE CORRECT//////
int check_initial_data(void)
{
  double xmax,aux_dens=0.0;
  
  xmax = 2.0*( 1.0+sqrt(1.0 + 3.0*(1.0+gamma*gamma)) )/(lambda*(1.0+gamma*gamma));
  
  if (p1<=0 || p2<=0 || p3<=0){
    
    printf("\n");
    printf("\t p1,p2 and p3 need to be positive. \n");
    printf("\n");
    exit(1);
  }
  
  if ( sqrt(p2*p3/(p1*p1*p1)) > xmax){
    printf("\n");
    printf("\t Wrong initial data for p1,p2,p3. \n");
    printf("\t The quantity sqrt(p2*p3/(p1*p1*p1)) need to be smaller than %lf. \n",xmax);
    printf("\n");
    exit(1);
  }

  if ( mu1c1 < -M_PI/2.0 || mu2c2 < -M_PI/2.0 || mu3c3 < -M_PI/2.0 ){
    
    printf("\n");
    printf("\t mu1c1, mu2c2, mu3c3 need to be greater than -pi/2 \n");
    printf("\n");
    exit(1);
  }

  if ( mu1c1 > 3.0*M_PI/2.0 || mu2c2 > 3.0*M_PI/2.0 || mu3c3 > 3.0*M_PI/2.0 ){
    
    printf("\n");
    printf("\t mu1c1, mu2c2, mu3c3 need to be smaller than 3pi/2 \n");
    printf("\n");
    exit(1);
  }
  
  if (tol_err < 0.0){
    printf("\n");
    printf("\t The tolerance error need to be positive. \n");
    printf("\n");
    exit(1);
  }
  
  if (tol_err < 1.0e-14){
    printf("\n");
    printf("\t The tolerance error need to be bigger than the machine error. \n");
    printf("\n");
    exit(1);
  }
  
  if (itmax < 0){
    printf("\n");
    printf("\t The maximal number of iterations need to be positive (and integer). \n");
    printf("\n");
    exit(1);
  }

  //SQRT ARGUMENT FOR THE MOMENTUM
  if (dyn_eq==0) {       //Classical    
    
    aux_dens = 1.0/(8.0*M_PI*G*gamma*gamma)
                 *(p1*p2*c1*c2 + p2*p3*c2*c3 + p1*p3*c1*c3 
		   + alpha*epsilon*p2*p3*c1 
		   - (1.0+gamma*gamma)*(alpha*p2*p3/(2.0*p1))
		   *(alpha*p2*p3/(2.0*p1)) 
		  ) -p1*p2*p3*potential();
    
  }else if (dyn_eq==1) { //Effective
        
    aux_dens = p1*p2*p3/(8.0*M_PI*G*gamma*gamma*Delta) 
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

  if (aux_dens < 0.0){
    printf("\n");
    printf("\t WARNING, the initial density is negative = %lf \n",aux_dens);
    printf("\n");
    exit(1);
  }

  return 0;
}


         ////////////////////
         /// IO ROUTINES/////
         ////////////////////

////////READ INITIAL CONDITIONS/////
int read_param(void)
{
  int log_pf;
  FILE *pf;

  pf=fopen(initfile,"r");
  
  log_pf=fscanf(pf,"%lf%lf%lf%lf%lf%lf %lf %lf%d %lf%lf%lf %d %d%lf%d  %d%d %d%d %lf%lf %lf%lf%lf",
		&mu1c1,&mu2c2,&mu3c3,&p1,&p2,&p3,
		&phi,
		&tol_err,&itmax,
		&initial_time,&final_time,&dt,
		&time_output,
		&Integrator,&alpha,&preferred_time,
		&dyn_eq,&std_out,
		&pot_switch, &pot_select,
		&m_phi, &l_phi,
		&V0_pot, &sigma1_pot, &sigma2_pot
		);

  fclose(pf);

  if (log_pf != 25) {
    printf("\n\t Error reading the parameters file.");
    printf("\n\t The program read %d initial values and they need to be 25. \n",log_pf);
    exit(1);  
  }

  return 0;
}


//// WRITE THE OBSERVABLES TO A FILE ////
int write_output(void)
{
  int log_pf;
  char newfile[200],*aux;  
  FILE *pf;

  if(density   < 1.0e-15)  density   = 0.0;
  if(shear     < 1.0e-15)  shear     = 0.0;
  if(Ricci     < 1.0e-15)  Ricci     = 0.0;
  if(fabs(expansion)  < 1.0e-15)  expansion = 0.0;
  if(fabs(constraint) < 1.0e-15) constraint = 0.0;
  
  aux="/density.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,density/density_crit);
  fclose(pf);

  aux="/shear.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,shear);
  fclose(pf);

  aux="/expansion.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,expansion);
  fclose(pf);

  aux="/Ricci.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,Ricci);
  fclose(pf);

  aux="/volume.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,volume);
  fclose(pf);

  aux="/a123.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \t %.15e \n",run_time,a_prom,a1,a2,a3);
  fclose(pf);

  aux="/H123.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \n",run_time,H1,H2,H3);
  fclose(pf);

  aux="/kasner.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \n",run_time,k1,k2,k3);
  fclose(pf);

  aux="/constraint.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \n",run_time,constraint);
  fclose(pf);

  aux="/constants.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e \t %.15e \t %.15e \n"
		 ,run_time,field_momentum(),omega,sigma2,curvature_param,c3*p3-c2*p2);
  fclose(pf);

  aux="/cp.t";
  strcpy(newfile,outputfile);
  strcat(newfile,aux);
  pf=fopen(newfile,"a");
  log_pf=fprintf(pf,"%.15e \t%.15e\t%.15e\t%.15e \t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",
		 run_time,c1,c2,c3,p1,p2,p3,mu1*c1,mu2*c2,mu3*c3);
  fclose(pf);

  if(pot_switch==1)
    {
      sphi = P_phi;
      if(preferred_time==1) sphi = sphi/volume;
	  	
      aux="/potential.t";
      strcpy(newfile,outputfile);
      strcat(newfile,aux);
      pf=fopen(newfile,"a");
      log_pf=fprintf(pf,"%.15e \t %.15e \t %.15e \t %.15e  \t %.15e \n"
		     ,run_time,phi,sphi,potential(),dev_pot());
      fclose(pf);
    }

  return log_pf=0;
}
