#include <iostream>
#include <fstream>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "Model_constants.h"

/**
 * \brief "computelib.h" defines main variables and all the functions used
 */

#define P2(x) ((x)*(x))
#define P3(x) ((x)*(x)*(x))


double *par; ///< par - model parameters
double *y; ///< y - model state variables
double *y_temp; ///< y_temp - temporary storage of y
double *y_converge; ///< y_converge - temporary storage of converged results

int *num_var; ///< num_var - total number of state variables
int *num_par; ///< num_par - total number of model parameters in macros file
int *num_cell; ///< num_cell - total number of cells coupled


double *L_cell; ///< Agonist concentration in the cells
double *L_rest; ///< Agonist cocentration before the stimulation starts

double *Q_ip3r; ///< Maximum rate of ip3r channel
double *Q_serca; ///< Maximum rate of serca pump
double *Q_ryr; ///< Maximum rate of ryr channel

void initialize(); ///< This function initializes the model variables and parameters
void allocatememory(int var_tot, int par_tot); ///< This function allocates memory for all the variables defined above
void file_names(char sub_folder, char timeseries_folder); ///< Creating FILE names of the output files
void singlecell(double tnow, double interval); ///< This function solves single cell dynamics
void dump_data(FILE *ft, FILE *fvm, FILE *fca,FILE *fpmca, FILE *fsrleak, FILE *fip3, FILE *fserca, FILE *fncx, FILE *fvocc,  FILE *fbkca, FILE *fip3r, FILE *fcacc,  FILE *fnsc, FILE *fsrca, FILE *fryr, FILE *fdag, double tnow); ///< This function writes output as .txt files
