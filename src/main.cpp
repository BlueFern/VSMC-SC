/**
 * \mainpage Vascular Smooth Muscle Cell Model
 * This model contains three fluid compartments, extracellular space (e), cytosol (cyt), and SR, which are separated by cell membrane and SR membrane.
 * Ion transport across cell membrane is given by four ion channels, large conductance calcium-activated potassium channel (BKCa),
 * voltage-operated calcium channel (VOCC), non-selective cation channel (NSC) and calcium-activated chloride channel (CACC),
 * one plasma membrane calcium ATPase pump (PMCA), and one Na/Ca exchanger (NCX).
 * Three calcium channels such as IP3 receptor channel (IP3R), Ryanodine receptor channel (RyR) and SR leak channel (SRleak),
 * and one SR calcium ATPase pump (SERCA) are included in the SR membrane. Calcium buffering in the cytosol is considered as a sum of calmodulin and other buffers.
 * Action of an agonist on adrenergic receptor cascade and the resulting formation of IP3 and DAG are included in the model.
 * @author Jaijus
 */


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>

#include "include/computelib.h"

int main(int argc, char* argv[]){

	/* Agonist concentration */
		double L_max = 301; // Maimum agonist concentration, nM
		double L_incre = 50; // increment of agonist concentration, nM
		double L_ini = 300.0; // initial agonist concentration, nM
		double L_var = L_ini; // Agonist concentration, nM

	/* Simulation controlling parameters */
		double tnow = 0.0; // Initial time , s
		double tfinal = 300; // Final time, s
		double interval = 1e-4; // Time interval, s
		double tol_state = 1e-4; // Tolerance limit for the convergence in the iteration
		double delay = 50; // Time allowed for variables to converge
		double stim_time = tfinal; // Total stimulation time, s
		double count_delay = 0; // control variable for the loop
		int tol_count = 0; // Control varialbe for fixed point iteration
		double file_write_freq = 0.01; // File writing frequency, s
		int folder_status; // variable to check the successful creartion of folder


		char main_folder[300];
		char sub_folder[300];
		char temp_name[300];

		/* Variables for coupled cells */
		num_var = (int *) malloc(sizeof(int)); // Total number of varibles
		num_par = (int *) malloc(sizeof(int)); // Total number of parameters
		num_cell = (int *) malloc(sizeof(int)); // Total number of cells

		num_var[0] 	= 18; // Total number of state variables in a cell
		num_par[0] 	= 40; // Total number of parameters in a cell

		int var_tot; // total number state variables in all the cells
		int par_tot; // total number of parameters in all the cells

		num_cell[0]	= 1; // Total number of SMCs

		var_tot	= num_cell[0] * num_var[0];
	 	par_tot	= num_cell[0] * num_par[0];

	 	allocatememory(var_tot, par_tot); // Allocating memory

		L_rest[0] = 0; // initial agonist concentration for the first 'delay' seconds


		/* Creating main folder */
		sprintf(main_folder,"%3.1fs_%1.4fs",tfinal,interval);
		folder_status = mkdir(main_folder,0777);
			if(folder_status == -1)
				{
				perror("Couldn't create sub-directory");
				return -1;
				}


	/* Loop to continue the simulation for the range of agonist concentration from L_ini to L_max   */
		while ( L_var < L_max)
		{
			/* Creating subfolder to save all the files for each loop */
			sprintf(sub_folder,"%s/full_%3.1fnM",main_folder,L_var);
			printf("%s\n",sub_folder);
			folder_status = mkdir(sub_folder,0777);
				if(folder_status == -1)
					{
					perror("Couldn't create sub-directory");
					return -1;
					}


			int count = 1; // counter to produce output files
			// OIpening and creating FILE names for output files

			FILE *fca; // Time series of intracellular calcium
			FILE *fvm; // Time series of membrane potential
			FILE *ft; // Time data
			FILE *fbkca; // Time series of BKCa channel flux
			FILE *fvocc; // Time series of BKCa channel flux
			FILE *fpmca; // Time series of PMCA channel flux
			FILE *fdag; // Time series of DAG concentration
			FILE *fserca; // Time series of SERCA pump flux
			FILE *fip3; // Time series of IP3 concentration
			FILE *fip3r; // Time series of IP3R channel flux
			FILE *fsrca; // Time series of intraSR calcium
			FILE *fsrleak; // Time series of SRleak flux
			FILE *fcacc; // Time series of CACC flux
			FILE *fnsc; // Time series of NSC flux
			FILE *fncx; // Time series of NCX channel flux
			FILE *fryr; // Time series of RyR channel flux

			sprintf(temp_name,"%s/Ca_timeseries.txt",sub_folder);
			fca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/Vm_timeseries.txt",sub_folder);
			fvm = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/time.txt",sub_folder);
			ft = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_BKCa.txt",sub_folder);
			fbkca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_VOCC.txt",sub_folder);
			fvocc = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_pmca.txt",sub_folder);
			fpmca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/DAG_timeseries.txt",sub_folder);
			fdag = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_serca.txt",sub_folder);
			fserca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/IP3_timeseries.txt",sub_folder);
			fip3 = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_IP3R.txt",sub_folder);
			fip3r = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/SR_Ca_timeseries.txt",sub_folder);
			fsrca = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_SRleak.txt",sub_folder);
			fsrleak = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_CACC.txt",sub_folder);
			fcacc = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_NSC.txt",sub_folder);
			fnsc = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_NCX.txt",sub_folder);
			fncx = fopen(temp_name,"w+");
			sprintf(temp_name,"%s/I_RyR.txt",sub_folder);
			fryr = fopen(temp_name,"w+");


			initialize(); // Initialize the variables
			dump_data(ft, fvm, fca, fpmca, fsrleak, fip3, fserca, fncx, fvocc, fbkca, fip3r, fcacc, fnsc , fsrca, fryr,fdag, tnow); // Creating output files

			// Simulation starts here
			while (tnow <= tfinal)
			{


				// Assign required agonist concentration to the cell after "delay" seconds.
				if (tnow >= (delay) && count_delay < 1)
				{
					count_delay = 2;
					for (int i = 0; i < num_cell[0]; i++)
						{
						L_cell[i] = L_var;
						}

				}


				while (tol_count<1)
				{
					tol_count = 2;
					singlecell(tnow, interval); // Solving single cell model equations

				/* Fixed point iteration checking*/
					for (int i=0;i<num_cell[0];i++)
					for (int j=0;j<num_var[0];j++)
						if (fabs(y_temp[i * num_var[0] + j] - y_converge[i * num_var[0] + j])>=tol_state)
						{
							tol_count = 0;
							/* Not converged */
							for (int k=0;k<num_cell[0];k++)
								for (int m=0;m<num_var[0];m++)
										y_converge[k * num_var[0] + m] = y_temp[k * num_var[0] + m];
							break;
						}
				}

				/* Converged*/
				/* Assigning converged results as the values at the new time step */
				for (int i=0;i<num_cell[0];i++)
					for (int j=0;j<num_var[0];j++)
						{
							y[i * num_var[0] + j] = y_temp[i * num_var[0] + j];
							y_converge[i * num_var[0] + j] = y_temp[i * num_var[0] + j];
						}

				tol_count = 0;

				tnow += interval; // new time step

				//  Checking file writing condition
					if (int(tnow/file_write_freq) == count)
						{
						dump_data(ft, fvm, fca, fpmca, fsrleak, fip3, fserca, fncx, fvocc, fbkca, fip3r, fcacc, fnsc , fsrca, fryr, fdag, tnow);
						count++;
						}

			}

			/* Closing all the files opened for writing */
			fclose(ft);
			fclose(fvm);
			fclose(fca);
			fclose(fpmca);
			fclose(fsrleak);
			fclose(fip3);
			fclose(fsrca);
			fclose(fncx);
			fclose(fvocc);
			fclose(fbkca);
			fclose(fserca);
			fclose(fip3r);
			fclose(fnsc);
			fclose(fcacc);
			fclose(fdag);
			fclose(fryr);
			tnow = 0.0;
			count_delay = 0;
			L_var += L_incre;
		}

	// free up memory
			free(y);
			free(y_temp);
			free(y_converge);
			free(par);
			free(num_var);
			free(num_par);
			free(num_cell);
			free(L_cell);
			free(L_rest);
			free(Q_ip3r);
			free(Q_serca);
			free(Q_ryr);

}
