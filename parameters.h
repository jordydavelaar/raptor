/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;

// METRIC PARAMETERS
////////////////////

// coordinate and metric choices
#define CAR  (0) // Minkowski
#define BL   (1) // Boyer-Lindquist,               x1=r, x2=th, and x3=phi
#define MBL  (2) // modified Boyer-Lindquist, x1=log(r), x2=th, and x3=phi
#define KS   (3) // Kerr-Schild,                   x1=r, x2=th, and x3=phi
#define MKS  (4) // modified Kerr-Schild,     x1=log(r), x2=th, and x3=phi
#define MKS2 (9) // Proper MKS coords

extern double a;
extern double R0; // Parameter for MKS coords

// Metric
#define metric (MKS2)
#if (metric == BL || metric == KS)
#define logscale (0) // Standard BL/KS coordinates; no logarithmic radius
#elif (metric == MBL || metric == MKS || metric == MKS2)
#define logscale (1) // Modified BL/KS coordinates; logarithmic radius
#endif

// MODEL PARAMETERS
///////////////////

// GRMHD data file
char GRMHD_FILE[256];
char OUTPUT_FILE[256];
int  SPHERICAL_ACC;
char TEMP_MODEL[100];

int ABSORPTION;
#define LIGHT_TRANSPORT                                                        \
    (1) // Toggle light transport calculation on/off for integration debugging

#define RT_OUTER_CUTOFF                                                        \
    (40.) // Outer boundary of radiative transfer computation

// Black hole mass
double  MBH;
double  M_UNIT;

// OBSERVER PARAMETERS
//////////////////////

double CAM_FREQ;
double TIME_INIT;
double INCLINATION;

// SED parameters
int    FREQS_PER_DEC;
double FREQ_MIN;
double FREQ_MAX;

#define source_dist (5.061e25) // Distance to M87 (cm); for Sgr A* use (2.47e22)
#define rcam (1e4) //(500.)    // Camera distance from the sing.(units of Rg)

int    IMG_WIDTH;
int    IMG_HEIGHT;
double CAM_SIZE_X;
double CAM_SIZE_Y;
#define max_order                                                              \
    (100) // Maximimum order of lensed image computed (0 = direct only)

// INTEGRATOR PARAMETERS
////////////////////////

#define OUTER_BOUND_POL (1000.) // Stop polarized integration beyond this radius

#define delta_num (1.e-6) // Used for numerical derivatives
#define max_steps (1e5)   // Maximum number of integration steps

double STEPSIZE;
#define cutoff_outer (1.1 * rcam) // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1.e-5) // Stop tracing at this distance from E.H. [BL]
#define VER (1)              //
#define RK4 (2)              //
#define int_method (2)       // method of integration 2=Verlet, 4-RK4

// MACROS
/////////

#define DIM 4
#define LOOP_i for (i = 0; i < DIM; i++)
#define LOOP_ij                                                                \
    for (i = 0; i < DIM; i++)                                                  \
        for (j = 0; j < DIM; j++)
#define LOOP_kl                                                                \
    for (k = 0; k < DIM; k++)                                                  \
        for (l = 0; l < DIM; l++)
#define LOOP_ijk                                                               \
    for (i = 0; i < DIM; i++)                                                  \
        for (j = 0; j < DIM; j++)                                              \
            for (k = 0; k < DIM; k++)
#define LOOP_ijkl                                                              \
    for (i = 0; i < DIM; i++)                                                  \
        for (j = 0; j < DIM; j++)                                              \
            for (k = 0; k < DIM; k++)                                          \
                for (l = 0; l < DIM; l++)

#endif // PARAMETERS_H
