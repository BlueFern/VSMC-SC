#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sys/stat.h>

/**
 * @file Memory_allocation.cpp
 * \brief This file contains "void allocatememory(int var_tot, int par_tot)" function
 */

extern double *y, *par, *y_temp, *y_converge;
extern int *num_cell;
extern double *Q_ip3r, *Q_serca, *Q_ryr;
extern double *L_cell, *L_rest;


/**************************************************************************************/
/**/void allocatememory(int var_tot, int par_tot) /**/
/**************************************************************************************/
{
	y = (double *) malloc(var_tot * sizeof(double));
	if (y == NULL) {
		printf("allocation failed at y\n");
		abort();
	}
	y_temp = (double *) malloc(var_tot * sizeof(double));
	if (y_temp == NULL) {
		printf("allocation failed at y_temp\n");
		abort();
	}
	y_converge = (double *) malloc(var_tot * sizeof(double));
	if (y_converge == NULL) {
		printf("allocation failed at y_converge\n");
		abort();
	}
	par = (double *) malloc(par_tot * sizeof(double));
	if (par == NULL) {
		printf("allocation failed at par\n");
		abort();
	}

	Q_ip3r = (double *) malloc(num_cell[0] * sizeof(double));
	if (Q_ip3r == NULL) {
		printf("allocation failed at Q_ip3r\n");
		abort();
	}

	Q_serca = (double *) malloc(num_cell[0] * sizeof(double));
	if (Q_serca == NULL) {
		printf("allocation failed at Q_serca\n");
		abort();
	}

	Q_ryr = (double *) malloc(num_cell[0] * sizeof(double));
	if (Q_ryr == NULL) {
		printf("allocation failed at Q_ryr\n");
		abort();
	}



	L_cell = (double *) malloc(num_cell[0] * sizeof(double));
	if (L_cell == NULL) {
		printf("allocation failed at L_cell\n");
		abort();
	}

	L_rest = (double *) malloc(sizeof(double));
}
