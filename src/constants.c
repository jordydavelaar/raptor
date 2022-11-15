/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////

double R_GRAV;        // Gravitational radius
double C_CONST;       // A constant frequently used in RTE integration
double R_ISCO;        // Innermost stable circular orbit
double CUTOFF_INNER;  // Inner integration boundary
double JANSKY_FACTOR; // Factor to scale image output

// FUNCTIONS
////////////

void set_constants() {
    // Horizon radius for integration cutoff
    double Rh = (1. + sqrt(1. - a * a));
    CUTOFF_INNER = Rh * (1. + horizon_marg); // Cutoff outside or inside BH EH
    R_GRAV = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm
    C_CONST = R_GRAV * PLANCK_CONSTANT /
              (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    // Innermost stable circular orbit (ISCO)
    double Z1 = 1. + pow(1. - a * a, (1. / 3.)) *
                         (pow(1. + a, (1. / 3.)) + pow(1. - a, (1. / 3.)));
    double Z2 = pow(3. * a * a + Z1 * Z1, 0.5);
    double retro =
        a < 0. ? 1. : -1.; // 1 for retrograde orbits, -1 for prograde
    R_ISCO = (3. + Z2 + retro * pow((3. - Z1) * (3. + Z1 + 2. * Z2), 0.5));

    JANSKY_FACTOR =
        1.e23; // 1.e23 is conversion from Jansky to ergs/Sr Hz s cm2
}
