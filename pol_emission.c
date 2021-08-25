/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "constants.h"
#include "functions.h"
#include "gsl/gsl_sf_hyperg.h"
#include "parameters.h"
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <stdio.h>

// POLARIZED COEFFICIENTS
/////////////////////////

// Dexter (2016) A.18
double I_I(double x) {
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) + 0.9977 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.19
double I_Q(double x) {
    return 2.5651 *
           (1 + 0.93193 * pow(x, -1. / 3.) + 0.499873 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.20
double I_V(double x) {
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
            0.0292545 * pow(x, -0.5) + 2.03773 * pow(x, -1. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.12
double j_I(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_I(x);
}

// A.13
double j_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_Q(x);
}

// A.14 (theta_B = pitch angle, k dot B)
double j_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / tan(theta_B) /
           3. / sqrt(3.) / SPEED_OF_LIGHT / theta_e / theta_e / theta_e *
           I_V(x);
}

// ABSORPTION
/////////////

// Use Kirchoff: a_nu = j_nu/B_nu

// ROTATION
///////////

// B.13 NOTE: make sure that, in B.13, Dexter does not mean X when he writes x
// (in the final term) Note that X = ... (B.8) UPDATE: it was apparently indeed
// a typo. while x = nu / nu_c UPDATE: it
// should be ln(X/120), not "ln x / 120". Thus I have greyed out the original,
// new one is below.
double f_m(double X) {
    return 2.011 * exp(-pow(X, 1.035) / 4.7) -
           cos(X * 0.5) * exp(-pow(X, 1.2) / 2.73) - 0.011 * exp(-X / 47.2) +
           (0.011 * exp(-X / 47.2) - pow(2., -1. / 3.) / pow(3., 23. / 6.) *
                                         10000. * M_PI * pow(X, -8. / 3.)) *
               0.5 * (1. + tanh(10. * log(X / 120.)));
}

// Bessel function approximations:
double bessel_appr(int n, double x) {
    if (n == 0)
        return -log(x / 2.) - 0.5772;

    if (n == 1)
        return 1. / x;

    if (n == 2)
        return 2. / x / x;

    exit(0);
}

// B.15:
double DeltaJ_5(double X) {
    return 0.43793091 * log(1. + 0.00185777 * pow(X, 1.50316886));
}

// B.4
double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    return 2. * M_PI * nu / 2. / SPEED_OF_LIGHT * wp2 * omega0 * omega0 /
           pow(2. * M_PI * nu, 4.) * f_m(Xe) *
           (bessel_appr(1, Thetaer) / bessel_appr(2, Thetaer) +
            6. * theta_e) *
           sin(theta_B) * sin(theta_B);
}

// B.14
double rho_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    return 2.0 * M_PI * nu / SPEED_OF_LIGHT * wp2 * omega0 /
           pow(2. * M_PI * nu, 3) * (bessel_appr(0, Thetaer) - DeltaJ_5(Xe)) /
           bessel_appr(2, Thetaer) * cos(theta_B);
}