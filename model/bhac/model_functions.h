/*
 * model_functions.h
 *
 * Written by J. Davelaar 2022
 *
 */

#ifndef MODEL_FUNCTIONS_H
#define MODEL_FUNCTIONS_H

void set_units(double);

void init_bhac_amr_data(char *fname);

void init_storage();

double interp_scalar_2D(double ***var, int i, int j, int k, double coeff[4]);

void Xtoij(double *X, int *i, int *j, double *del);

double interp_scalar(double **var, int c, double coeff[4]);

void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);

double get_detgamma(double x, double y, double z);

int find_igrid(double x[4], struct block *block_info, double ***Xc);

#endif
