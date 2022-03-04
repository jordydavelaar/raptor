/*
 * raptor_harm_model.h
 *
 * Written by J. Davelaar
 *
 */

#define NDIM 4
#define NPRIM 8

struct GRMHD {
    double U_u[4];
    double B_u[4];
    double U_d[4];
    double B_d[4];
    double n_e;
    double B;
    double theta_e;
    double sigma;
    double dx_local;
    int igrid_c;
};

int KRHO = 0;
int UU = 1;
int U1 = 2;
int U2 = 3;
int U3 = 4;
int B1 = 5;
int B2 = 6;
int B3 = 7;

double ****p;
double R_HIGH, R_LOW;

double hslope = 1;

//		     0 1  2  3  4   5  6  7  8  9    10   11
// wnames          = 'd s1 s2 s3 tau b1 b2 b3 Ds dtr1 lfac xi'
int D = 0;
int S1 = 1;
int S2 = 2;
int S3 = 3;
int TAU = 4;
int DS = 8;
int LFAC = 10;
int XI = 11;

double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Ne_unit;
double Thetae_unit;

double R0;
double stopx[4], startx[4], dx[4];
double a;
double th_len;
int NP = 8;
int N1, N2, N3;

double *neqpar;
double xprobmin[3], xprobmax[3];
int ng[3], nxlone[3], nleafs;
double *eqpar;

struct block {
    int ind[3], level, size[3];
    double lb[3], dxc_block[3];
};

int count_leaf = 0, count_node = 0;
int *forest; //= NULL;
struct block *block_info;
int block_size = 1;
int forest_size = 1;

#ifndef RAPTOR_HARM_MODEL_H
#define RAPTOR_HARM_MODEL_H

#include "functions.h"
#include "parameters.h"
#include <stdio.h>

double Ladv, dMact;

double **Xcoord;
double ***Xgrid;
double ***Xbar;

// RAPTOR_HARM_MODEL.C
//////////////////////

// See grmonty paper by Dolence et al.
// HARM model internal utilities
void set_units(double);
void init_bhac_amr_data(char *fname);
void init_storage();
double interp_scalar_2D(double ***var, int i, int j, int k, double coeff[4]);
void Xtoij(double *X, int *i, int *j, double *del);

void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);

double get_detgamma(double x, double y, double z);

#endif // RAPTOR_HARM_MODEL_H
