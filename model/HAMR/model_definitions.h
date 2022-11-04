/*
 * model_definitions.h
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
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
