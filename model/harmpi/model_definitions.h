/*
 * model_definitions.h
 *
 * harmpi model for RAPTOR
 *
 * Based on public harmpi code: https://github.com/atchekho/harmpi
 *
 * Authors Anna Chashkina, Jordy Davelaar
 *
 */

#ifndef MODEL_DEFINITIONS_H
#define MODEL_DEFINITIONS_H

#define NDIM 4
#define NPRIM 8

typedef struct GRMHD {
    double U_u[4];
    double B_u[4];
    double U_d[4];
    double B_d[4];
    double n_e;
    double B;
    double sigma;
    double sigma_min;
    double beta;
    double theta_e;
    double dx_local;
    int igrid_c;
} GRMHD;

/* mnemonics for primitive vars; conserved vars */
#define KRHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7

#endif // MODEL_DEFINITIONS_H
