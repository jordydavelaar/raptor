/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "constants.h"
#include "parameters.h"
#include <math.h>
#include <stdio.h>

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

    // Calibration constant for the spectral irradiance
    // We want Jansky/pixel^2.
    double d_x = CAM_SIZE_X * R_GRAV; // Size of image in cm
    double delta_x =
        d_x / source_dist; // it is far away so just a ratio without tan
    double d_y = CAM_SIZE_Y * R_GRAV; // Size of image in cm
    double delta_y =
        d_y / source_dist; // it is far away so just a ratio without tan
    double pixsize = (delta_x / (double)IMG_WIDTH) *
                     (delta_y / (double)IMG_HEIGHT); // Pix size in Sr
    JANSKY_FACTOR =
        1.e23 * pixsize; // 1.e23 is conversion from Jansky to ergs/Sr Hz s cm2
}
