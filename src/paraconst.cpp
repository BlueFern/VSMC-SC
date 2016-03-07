#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>
#include "include/macros.h"
#include "include/Model_constants.h"

/**
 * @file paraconst.cpp
 * \brief This file contains "void paraconst(FILE *fparaconst, double tfinal, double interval, double file_write_freq, double delay, double stim_time, double tol_state)" function
 * The paraconst function write initial values of state variables, all the constant and paramters of the model to Readme.txt file.
 */

extern double *L_rest;
extern double *y, *par;
extern double *L_rest;
extern int *num_par, *num_var;

void paraconst(FILE *fparaconst, double tfinal, double interval, double file_write_freq, double delay, double stim_time, double tol_state)

	{
	fprintf(fparaconst,"Simulation controlling parameters\n" );
	fprintf(fparaconst,"Final time, t_final  = %4.1f s\n",  tfinal );
	fprintf(fparaconst,"Time interval, interval  = %f s\n",  interval );
	fprintf(fparaconst,"Tolerance limit, tol_state  = %f s\n",  tol_state );
	fprintf(fparaconst,"File writing frequency, file_write_freq  = %2.2f s-1\n",  file_write_freq );
	fprintf(fparaconst,"Time delay for converging variables and functions, delay  = %4.1f s\n", delay );
	fprintf(fparaconst,"Stimulation time, stim_time  = %4.1f s\n",  stim_time );

	fprintf(fparaconst, "\nSimulation constants\n" );
	fprintf(fparaconst,"Agonist concentration at rest, L_rest  = %3.1f nM\n",  L_rest[0] );

	fprintf(fparaconst, "\nStandard VSMC Parameters\n" );
	fprintf(fparaconst,"Cell capacitance, C_m  = %3.1f  pC mV-1\n",  C_m );
	fprintf(fparaconst,"Faraday's constant, F  = %5.1f  C mol-1\n",  F );
	fprintf(fparaconst,"Universal Gas constant, R  = %4.1f  mJ mol-1K-1\n",  R );
	fprintf(fparaconst,"Temperature, T  = %3.1f  K\n",  T );
	fprintf(fparaconst,"Avagadro's Number  = %e  \n",  N_AV );
	fprintf(fparaconst,"VSMC volume, Vol_smc  = %2.1f  pL\n",  Vol_smc );
	fprintf(fparaconst,"Cytosol volume, Vol_cyt  = %2.1f  pL\n",  Vol_cyt );
	fprintf(fparaconst,"Sarcoplasmic reticulum volume, Vol_SR  = %2.2f  pL\n",  Vol_SR );
	fprintf(fparaconst,"Ion concentration to current conversion factor for cytosol, alpha  = %e  pC nM-1\n",  alpha );
	fprintf(fparaconst,"Ion concentration to current conversion factor for intraSR, beta  = %e  pC nM-1\n",  beta );
	fprintf(fparaconst,"Potassium Concentration in cytosol, K_cyt  = %1.2e  nM\n",  K_cyt );
	fprintf(fparaconst,"Sodium Concentration in cytosol, Na_cyt  = %1.2e  nM\n",  Na_cyt );
	fprintf(fparaconst,"Chloride Concentration in cytosol, Cl_cyt  = %1.2e  nM\n",  Cl_cyt );
	fprintf(fparaconst,"Calcium Concentration in extracellular space, Ca_e  = %1.2e  nM\n",  Ca_e );
	fprintf(fparaconst,"Potassium Concentration in extracellular space, K_e  = %1.2e  nM\n",  K_e );
	fprintf(fparaconst,"Sodium Concentration in extracellular space, Na_e  = %1.2e  nM\n",  Na_e );
	fprintf(fparaconst,"Chloride Concentration in extracellular space, Cl_e  = %1.2e  nM\n",  Cl_e );
	fprintf(fparaconst,"Nernst potential of potassium, E_K  = %3.3f  mV\n",  par[ 0*num_par[0] + E_K] );
	fprintf(fparaconst,"Nernst potential of Sodium, E_Na  = %3.3f  mV\n",  par[ 0*num_par[0] + E_Na] );
	fprintf(fparaconst,"Nernst potential of chloride, E_Cl  = %3.3f  mV\n",  par[0*num_par[0] + E_Cl] );
	fprintf(fparaconst,"\nInitial values of state variables\n");
	fprintf(fparaconst,"Initial Concentration of intraSR calcium, SR_Ca_ini = %1.4e  nM\n",  SR_Ca_ini );
	fprintf(fparaconst,"Initial Concentration of intracellular calcium, Ca_ini = %3.2f nM\n", Ca_ini);
	fprintf(fparaconst,"Initial Concentration of Membrane potential, vm_ini  = %3.2f mV\n", vm_ini);

	fprintf(fparaconst,"\nParameters for Membrane Channels, Pumps and exchangers\n");
	fprintf(fparaconst,"Whole cell conductance of VOCC  = %1.3e  nS\n",  G_vocc );
	fprintf(fparaconst,"Whole cell conductance of BKCa channel  = %1.3e  nS\n",  G_bkca );
	fprintf(fparaconst,"Fast activation time constant  = %f  s\n",  t_pf );
	fprintf(fparaconst,"Slow activation time constant  = %f  s\n",  t_ps );
	fprintf(fparaconst,"Whole cell conductance of CACC  = %1.2f  nS\n",  G_cacc );
	fprintf(fparaconst,"order of the function shows depenedncy on intracellular calcium  = %1.1f  \n",  n_cacc );
	fprintf(fparaconst,"Dissociation constant of CACC  = %3.1f  nM\n",  k_cacc );
	fprintf(fparaconst,"Whole cell conductance of NCX channel  = %1.2e  pS\n",  G_ncx );
	fprintf(fparaconst,"Position of the energy barrier of NCX channel  = %1.2f  nM\n",  gamma_ncx );
	fprintf(fparaconst,"Dissociation constant for allosteric component - NCX channel  = %3.1f  nM\n",  k_ncx );
	fprintf(fparaconst,"d_ncx  = %1.3f  \n",  d_ncx );
	fprintf(fparaconst,"Reverse potential of NSC channel  = %1.1f  mV\n",  E_nsc );
	fprintf(fparaconst,"Dissociation constant of NSC  = %4.1f  nM\n",  k_nsc );
	fprintf(fparaconst,"Whole cell conductance of NSC  = %1.2f  nS\n",  G_nsc );
	fprintf(fparaconst,"Dissociation constant of PMCA  = %3.1f  nM\n",  k_pmca );
	fprintf(fparaconst,"Maximum rate of calcium flux of PMCA  = %1.4e  nM s-1\n",  Q_pmca );
	fprintf(fparaconst,"Order of calcium dependency, PMCA pump  = %1.1f  nM\n",  n_pmca );

	fprintf(fparaconst,"\nParameters for Sarcoplasmic reticulum channels and pumps\n");
	fprintf(fparaconst,"Binding Constants for IP3R channel\n");

	fprintf(fparaconst,"\nInitial values\n");
for (int i = 0; i < num_var[0] ; i++)
	fprintf(fparaconst,"y[%d]  = %f  \n", i, y[ 0*num_var[0] + i] );

	for (int i=0; i<num_par[0]; i++)
	{
	fprintf(fparaconst,"par[%d] = %1.10e\n", i, par[0 * num_par[0] + i]);
	}
	}
