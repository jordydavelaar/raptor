/*
 * raptor_harm_model.h
 *
 * Written by J. Davelaar
 *
 */

#include "functions.h"
#include "parameters.h"

#ifndef RAPTOR_HARM_MODEL_H
#define RAPTOR_HARM_MODEL_H

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

#define KRHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7
#define LF 8

double ****p;
double R_HIGH, R_LOW;

double hslope;

//		     0 1  2  3  4   5  6  7  8  9    10   11
// wnames          = 'd s1 s2 s3 tau b1 b2 b3 Ds dtr1 lfac xi'
#define D 0
#define S1 1
#define S2 2
#define S3 3
#define TAU 4
#define DS 8
#define LFAC 10
#define XI 11

double L_unit;
double T_unit;
double RHO_unit;
double U_unit;
double B_unit;
double Ne_unit;
double Thetae_unit;

double R0;
double stopx[4], startx[4], dx[4];
double a, Q;
double th_len;
#define NP 8;
int N1, N2, N3;

double *neqpar;
double xprobmin[3], xprobmax[3];
int ng[3], nxlone[3], nleafs;
double *eqpar;

struct block {
    int ind[3], level, size[3];
    double lb[3], dxc_block[3];
};

int *forest;
int *nx;
struct block *block_info;
int block_size;
int forest_size;
int cells;
double buffer[1];
unsigned int buffer_i[1];

double **values;
long int offset;

int levmaxini, ndimini, ndirini, nwini, nws, neqparini, it, t;

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
