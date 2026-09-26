/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <complex.h>
#include <hdf5.h>
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

#define IMGFILE (1)
#define SPECFILE (1)
#define RAD_TRANS (1)
#define POL (1)

// Polarized-transport and geodesic accuracy diagnostics in the image file
// (fnorm_*, wp_*, geo_*): 1 to enable (adds the Camera fields below)
#define POL_DIAGNOSTICS (0)

#define num_frequencies 1

#define FREQFILE (0)
#define FREQLOG (1)

#define FREQS (FREQLOG)

#define EMISUSER (0)

#define UNIF (1)

#define AMR 0
#define SMR 1

#define num_pixels_1d 10
#define tot_pixels 100

#define USERSPEC (1)
#define nspec 4

typedef struct Camera {
    double IQUV[tot_pixels][num_frequencies][4]; // intensity
    double tau[tot_pixels][num_frequencies];     // intensity
    double tauF[tot_pixels][num_frequencies];    // intensity
#if (POL_DIAGNOSTICS)
    double fnorm_dev[tot_pixels][num_frequencies]; // max |f_perp|^2 - 1 on ray
    double fnorm_th[tot_pixels][num_frequencies];  // theta at that max
    double fnorm_cam[tot_pixels][num_frequencies]; // |f_perp|^2 at camera
    double wp_dchi[tot_pixels][num_frequencies]; // Walker-Penrose EVPA error (deg)
    double wp_dpsi[tot_pixels][num_frequencies]; // Poincare-sphere error (deg)
#endif
    double alpha[tot_pixels];                    // impact parameter
    double beta[tot_pixels];                     // impact parameter
    double norder[tot_pixels]; // winding number n = Delta_phi / (2*pi)
    double norder_mino[tot_pixels]; // Mino-time half-orbit count (GL20/Himwich)
    int ncross[tot_pixels];         // number of equatorial-plane crossings
    double kappa1[tot_pixels];      // Walker-Penrose constant at camera (Re)
    double kappa2[tot_pixels];      // Walker-Penrose constant at camera (Im)
    double dchi_grav[tot_pixels];   // EVPA rotation, source at infinity -> camera
#if (POL_DIAGNOSTICS)
    double geo_dlam[tot_pixels];   // max |lambda - lambda_cam| along ray
    double geo_deta[tot_pixels];   // max |eta - eta_cam| along ray
    double geo_thpole[tot_pixels]; // closest approach to polar axis (rad)
#endif
    double lcorner[2];                           // lower left corner of a block
    double dx[2];                                // pixel spacing of block
    int level;
    int ind[2];
} Camera;

#define sign(x) (((x) < 0) ? -1 : ((x) > 0))

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

#define SIGMA_CUT (1)
#define THETAE_MAX (100)
#define THETAE_MIN (1e-3)

#define LIGHT_TRANSPORT                                                        \
    (1) // Toggle light transport calculation on/off for integration debugging

// OBSERVER PARAMETERS
//////////////////////
#define rcam (1e4) //(500.)    // Camera distance from the sing.(units of Rg)

#define max_order (100) // Maximimum order of lensed images 0 = direct only

// INTEGRATOR PARAMETERS
////////////////////////

#define RT_OUTER_CUTOFF (50.) // Stop polarized integration beyond this radius

#define delta_num (1.e-4) // Used for numerical derivatives
#define max_steps (1e5)   // Maximum number of integration steps

#define cutoff_outer (1.1 * rcam) // Outer cutoff, near flat spacetime, in M
#define horizon_marg (1.e-2) // Stop tracing at this distance from E.H. [BL]
#define RK2 (1)              //
#define VER (2)              //
#define RK4 (3)              //
#define RK45 (4)             //
#define int_method (RK4)

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

#endif // DEFINITIONS_H
