/*
 * model_definitions.h
 *
 * Written by J. Davelaar 2022
 *
 */

#ifndef MODEL_DEFINITIONS_H
#define MODEL_DEFINITIONS_H

#define SFC 1

#define KRHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7

#define D 0
#define S1 1
#define S2 2
#define S3 3
#define TAU 4
#define DS 8
extern int LFAC;
extern int XI;

#define NP 8

#define NSPIN 3

#define PLUS 1
#define MINUS -1
#define BPOL (PLUS)

typedef struct GRMHD {
    double U_u[4];
    double B_u[4];
    double U_d[4];
    double B_d[4];
    double n_e;
    double B;
    double theta_e;
    double sigma;
    double sigma_min;
    double dx_local;
    double beta;
    int igrid_c;
} GRMHD;

typedef struct block {
    int ind[3], level, size[3];
    double lb[3], dxc_block[3];
} block;

#endif
