/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 * RAPTOR uses cgs units for light transport calculations.
 * Entries marked [BL] are only applicable to Boyer-Lindquist coordinates.
 */

#ifndef GLOBAL_VARS_H
#define GLOBAL_VARS_H

// CAMERA.C
///////////

extern int num_blocks, tot_blocks, tot_blocks_alloc;
extern double BLOCK_SIZE_X, BLOCK_SIZE_Y;
extern int max_level;

// CORE.C
/////////

extern char GRMHD_FILE[256];

extern double MBH, M_UNIT, TIME_INIT, INCLINATION;
extern double R_HIGH, R_LOW;
extern double FREQS_PER_DEC, FREQ_MIN, FREQ_MAX;

extern double SOURCE_DIST; // Distance to M87 (cm); for Sgr A* use (2.47e22)

extern int IMG_WIDTH, IMG_HEIGHT;
extern double CAM_SIZE_X, CAM_SIZE_Y;
extern double STEPSIZE;

// CONSTANTS.C
//////////////

extern double R_GRAV;        // Gravitational radius
extern double C_CONST;       // A constant frequently used in RTE integration
extern double R_ISCO;        // Innermost stable circular orbit
extern double CUTOFF_INNER;  // Inner integration boundary
extern double JANSKY_FACTOR; // Factor to scale image output

#endif
