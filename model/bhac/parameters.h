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

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEBUG (1)

#define NDIM 4
#define NPRIM 8

char inputfile[100];
#define IMGFILE (1)
#define SPECFILE (1)
#define RAD_TRANS (1)
#define POL (0)

#define AMR 0
#define SMR 1

#define num_frequencies 1
#define num_pixels_1d 10
#define tot_pixels 100
int num_blocks, tot_blocks, tot_blocks_alloc;
double BLOCK_SIZE_X, BLOCK_SIZE_Y;
int max_level;

struct Camera {
    double IQUV[tot_pixels][num_frequencies][4]; // intensity
    double tau[tot_pixels][num_frequencies];     // intensity
    double alpha[tot_pixels];                    // impact parameter
    double beta[tot_pixels];                     // impact parameter
    double lcorner[2];                           // lower left corner of a block
    double dx[2];                                // pixel spacing of block
    int level;
    int ind[2];
};

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

extern double L_unit;
extern double T_unit;
extern double RHO_unit;
extern double U_unit;
extern double B_unit;
extern double Ne_unit;
extern double Thetae_unit;

// DISTRIBUTION CHOISES
//////////////////
#define KAPPA (0) // kappa distribution
#define TH (1)    // thermal distribution
#define POWER (2) // Power-law distribution
#define DF (TH)   // Distribution function
#define kappa 5.0
#define power 2.5
#define gamma_min 1.
#define gamma_max 1000.

// METRIC PARAMETERS
////////////////////

// coordinate and metric choices
#define CAR (0) // Minkowski
#define BL (1)  // Boyer-Lindquist,               x1=r, x2=th, and x3=phi
#define MBL (2) // modified Boyer-Lindquist, x1=log(r), x2=th, and x3=phi
#define KS (3)  // Kerr-Schild,                   x1=r, x2=th, and x3=phi
#define MKS (4) // modified Kerr-Schild,          x1=log(r), x2=th, and x3=phi
#define MKSHARM                                                                \
    (5) // HARM3D MKS coords        x1=log(r), x2=th/pi, and x3=phi/2pi
#define MKSBHAC                                                                \
    (6) // BHAC style MKS coords          x1=log(r), x2=th/pi, and x3=phi
#define MKSN (7) //  modified Kerr-Schild-Newman,  x1=log(r), x2=th, and x3=phi
#define CKS (8)  //  modified Kerr-Schild-Newman,  x1=log(r), x2=th, and x3=phi

extern double a, Q;
extern double R0; // Parameter for MKS coords

// Metric
#define metric (MKSBHAC)
#if (metric == BL || metric == KS || metric == CKS)
#define logscale (0) // Standard BL/KS coordinates; no logarithmic radius
#elif (metric == MBL || metric == MKS || metric == MKSHARM ||                  \
       metric == MKSBHAC || metric == MKSN)
#define logscale (1) // Modified BL/KS coordinates; logarithmic radius
#endif

// MODEL PARAMETERS
///////////////////

// GRMHD data file
char GRMHD_FILE[256];
char OUTPUT_FILE[256];
int SPHERICAL_ACC;
char TEMP_MODEL[100];

int ABSORPTION;
#define LIGHT_TRANSPORT                                                        \
    (1) // Toggle light transport calculation on/off for integration debugging

#define RT_OUTER_CUTOFF                                                        \
    (40.) // Outer boundary of radiative transfer computation

// Black hole mass
double MBH;
double M_UNIT;

// OBSERVER PARAMETERS
//////////////////////

double CAM_FREQ;
double TIME_INIT;
double INCLINATION;

// SED parameters
int FREQS_PER_DEC;
double FREQ_MIN;
double FREQ_MAX;

double source_dist; // Distance to M87 (cm); for Sgr A* use (2.47e22)
#define rcam (1e4)  //(500.)    // Camera distance from the sing.(units of Rg)

int IMG_WIDTH;
int IMG_HEIGHT;
double CAM_SIZE_X;
double CAM_SIZE_Y;
#define max_order                                                              \
    (100) // Maximimum order of lensed image computed (0 = direct only)

// INTEGRATOR PARAMETERS
////////////////////////

#define OUTER_BOUND_POL (1000.) // Stop polarized integration beyond this radius

#define delta_num (1.e-4) // Used for numerical derivatives
#define max_steps (1e5)   // Maximum number of integration steps

double STEPSIZE;
#define cutoff_outer (1.1 * rcam) // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1.e-5) // Stop tracing at this distance from E.H. [BL]
#define VER (1)              //
#define RK4 (2)              //
#define RK45 (3)             //
#define int_method (2)       // method of integration 2=Verlet, 4-RK4

// CONSTANTS
////////////

// PHYSICAL CONSTANTS
#define ELECTRON_CHARGE (4.80320425e-10)
#define ELECTRON_MASS (9.1093829e-28)
#define PROTON_MASS (1.6726219e-24)
#define BOLTZMANN_CONSTANT (1.3806488e-16)
#define SPEED_OF_LIGHT (2.99792458e10)
#define PLANCK_CONSTANT (6.62606885e-27)
#define MPCL2 (0.0015033)
#define GGRAV (6.674e-8)
#define MSUN (1.989e33)
#define KPCTOCM (3.086e21)
#define MPoME (PROTON_MASS / ELECTRON_MASS)
//#define M_PI 3.14159265358979323846

// Constants that must be evaluated at startup
// (They depend on spin and other user-supplied parameters)
double R_GRAV;        // Gravitational radius
double C_CONST;       // A constant frequently used in RTE integration
double R_ISCO;        // Innermost stable circular orbit
double CUTOFF_INNER;  // Inner integration boundary
double JANSKY_FACTOR; // Factor to scale image output

// MACROS
/////////

#define DIM 4
#define LOOP_i for (int i = 0; i < DIM; i++)
#define LOOP_ij                                                                \
    for (int i = 0; i < DIM; i++)                                              \
        for (int j = 0; j < DIM; j++)
#define LOOP_kl                                                                \
    for (int k = 0; k < DIM; k++)                                              \
        for (int l = 0; l < DIM; l++)
#define LOOP_ijk                                                               \
    for (int i = 0; i < DIM; i++)                                              \
        for (int j = 0; j < DIM; j++)                                          \
            for (int k = 0; k < DIM; k++)
#define LOOP_ijkl                                                              \
    for (int i = 0; i < DIM; i++)                                              \
        for (int j = 0; j < DIM; j++)                                          \
            for (int k = 0; k < DIM; k++)                                      \
                for (int l = 0; l < DIM; l++)

#endif // PARAMETERS_H
