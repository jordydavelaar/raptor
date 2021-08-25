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

// kappa distribution function

int N_theta, N_theta_e, N_B, N_n_e, N_nuratio, N_nu;
double theta_min, theta_max, d_theta;
double theta_e_min, theta_e_max, d_theta_e;
double B_min, B_max, d_B;
double n_e_min, n_e_max, d_n_e;
double nuratio_min, nuratio_max, d_nuratio;
double nu_min, nu_max, d_nu;

double kappa, gamma_min, gamma_max, gamma_cutoff;
double **j_nu_data, **alpha_nu_data;

#define ME (9.10956e-28)
#define mc2 (8.18726e-07)
#define kb (1.38e-16)
#define hpl (6.6262e-27)
#define CL (2.99792458e10)
#define keV (1.602e-9)
#define alphaf (7.29735e-3)
#define h__mc2 (8.09e-21)
#define SIGMATH (0.665245873e-24)

double Xmax = 1e-25, Xmin = 1e25;

// non thermal emission
double emission_coeff_kappa_FIT(double nu, double Ne, double Thetae, double B,
                                double theta) {
    // emissivity for the kappa distribution function, see Pandya et al. 2016
    double nuc, sth, nus, x, w, X_kappa, factor;
    double J_low, J_high, J_s;

    double Rhigh = 4.5;
    double Rlow = 4.5;
    double b2 = pow((B / B_unit) / (Ne / Ne_unit), 2.);
    // printf("%g\n", b2);
    kappa = Rhigh * b2 / (1 + b2) + Rlow / (1 + b2);

    w = Thetae; // sqrt(  2./9./kappa *Thetae * Thetae);
    nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    sth = sin(theta);

    factor = (Ne * pow(ELECTRON_CHARGE, 2.) * nuc * sth) / SPEED_OF_LIGHT;

    //      fprintf(stderr,"sinth %g\n", sth);
    nus = nuc * sth * pow(w * kappa, 2);
    if (nu > 1.e12 * nus || Thetae > 400. || Thetae < .1)
        return (0.);
    X_kappa = nu / nus;
    //      fprintf(stderr, "X_kappa %g\n", X_kappa);
    J_low = pow(X_kappa, 1. / 3.) * sth * 4 * M_PI * tgamma(kappa - 4. / 3.) /
            (pow(3, 7. / 3.) * tgamma(kappa - 2.));
    //      fprintf(stderr, "J_low %g\n", J_low);
    J_high = pow(X_kappa, -(kappa - 2) / 2.) * sth * pow(3, (kappa - 1) / 2.) *
             (kappa - 2.) * (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
             tgamma(kappa / 4. + 4. / 3.);
    //      fprintf(stderr, "J_high %g\n", J_high );
    x = 3 * pow(kappa, -3. / 2.);

    J_s = pow((pow(J_low, -x) + pow(J_high, -x)), -1. / x);
    //      fprintf(stderr,"J_s %g\n", J_s * factor);
    return (J_s * factor);
}

double absorption_coeff_kappa_FIT(double nu, double Ne, double Thetae, double B,
                                  double theta) {
    // absortion for the kappa distribution function, see Pandya et al. 2016
    double nuc, sth, nus, x, w, X_kappa, factor;
    double A_low, A_high, A_s;
    double Rhigh = 4.5;
    double Rlow = 4.5;
    double b2 =
        pow((B / B_unit) / (Ne / Ne_unit), 2.); //(B/B_unit)/(Ne/Ne_unit);
    kappa = Rhigh * b2 / (1 + b2) + Rlow / (1 + b2);
    w = Thetae; // sqrt(2./9./kappa *Thetae * Thetae);
    nuc = ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    sth = sin(theta);

    factor = Ne * ELECTRON_CHARGE / (B * sth);

    nus = nuc * sth * pow(w * kappa, 2);
    if (nu > 1.e12 * nus || Thetae > 400. || Thetae < 1.)
        return (0.);
    X_kappa = nu / nus;

    // identity to be able to calculate a hypergeometric function, from the code
    // symphony by Pandya et al. 2016
    double a = kappa - 1. / 3.;
    double b = kappa + 1.;
    double c = kappa + 2. / 3.;
    double z = -kappa * w;
    double hyp2F1 = pow(1. - z, -a) * tgamma(c) * tgamma(b - a) /
                        (tgamma(b) * tgamma(c - a)) *
                        gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z)) +
                    pow(1. - z, -b) * tgamma(c) * tgamma(a - b) /
                        (tgamma(a) * tgamma(c - b)) *
                        gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1. - z));

    A_low = pow(X_kappa, -5. / 3.) * pow(3, 1. / 6.) * (10. / 41.) *
            pow(2 * M_PI, 2) / pow(w * kappa, 16. / 3. - kappa) * (kappa - 2.) *
            (kappa - 1.) * kappa / (3. * kappa - 1.) * tgamma(5. / 3.) * hyp2F1;
    //      fprintf(stderr, "A_low %g\n", A_low);
    A_high = pow(X_kappa, -(3. + kappa) / 2.) * (pow(M_PI, 5. / 2.) / 3.) *
             ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 5.)) *
             (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) *
             (pow(3. / kappa, 19. / 4.) + 3. / 5.);

    //      fprintf(stderr, "A_high %g\n", A_high);

    x = pow(-7. / 4. + 8. * kappa / 5., -43. / 50.);

    //      fprintf(stderr,"factor %g\n",factor);

    A_s = pow((pow(A_low, -x) + pow(A_high, -x)), -1. / x);
    //      fprintf(stderr, "A_s %g\n", A_s*factor);
    return (factor * A_s);
}


// Return emissivity j_nu which depends on local plasma parameters
// Ref. Dolence & Moscibrodzka 2009
double emission_coeff_THSYNCH(double B, double theta, double THETA_e,
                              double nu_plasma, double n_e) {
    double nu_c =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);

    double nu_s = 2. / 9. * nu_c * THETA_e * THETA_e * sin(theta);

    double X = nu_plasma / (nu_s);
    double f = pow(pow(X, 0.5) + pow(2., 11. / 12.) * pow(X, 1. / 6.), 2.);
    double j_nu = n_e * sqrt(2.) * M_PI * ELECTRON_CHARGE * ELECTRON_CHARGE *
                  nu_s / (6. * THETA_e * THETA_e * SPEED_OF_LIGHT) * f *
                  exp(-pow(X, 1. / 3.));

    return j_nu;
}

// Return emission constant j_nu as described in Dexter (2009) (the geokerr
// paper)
double emission_coeff_THSYNCHAV(double B, double THETA_e, double nu_p,
                                double n_e) {
    double nu_c = 3. * ELECTRON_CHARGE * B * THETA_e * THETA_e /
                  (4. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double x_M = nu_p / nu_c;
    double I = 4.0505 / pow(x_M, 1. / 6.) *
               (1. + 0.4 / pow(x_M, 1. / 4.) + 0.5316 / sqrt(x_M)) *
               exp(-1.8899 * pow(x_M, 1. / 3.));

    double j_nu = nu_p * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                  (2. * sqrt(3.) * SPEED_OF_LIGHT) * 1. / (THETA_e * THETA_e) *
                  I;

    //    j_nu = THETA_e * (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
    //    BOLTZMANN_CONSTANT; j_nu = n_e;

    return j_nu;
}

// Return emissivity for the simple Gaussian hot spot model discussed in Dexter
// 2009.
double emissivity_hotspot(double *X_u) {
    double xspot[4];

    double r = logscale ? exp(X_u[1]) : X_u[1];

    double Rspot = 0.5; // Spot size (radius)

    double r_spot = 6.0; // Spot orbit radius
    double th_spot = 0.5 * M_PI;
    double r32 = pow(r_spot, 1.5);
    double omega = 1. / (r32 + a);
    double P = 2. * M_PI / omega; // Period of the spot on a Keplerian orbit[M]

    // spot currrent position
    xspot[0] = X_u[0]; // current coordinate time
    xspot[1] = logscale ? log(r_spot) : r_spot;
    xspot[2] = th_spot; // equator 0.5*pi
    xspot[3] =
        fmod(X_u[0] / P, 1.) * 2. * M_PI + M_PI; // spot current phi at t=X[0]

    // Pseudo-Cartesian coordinates
    double xc = sqrt(r * r + a * a) * cos(X_u[3]);
    double yc = sqrt(r * r + a * a) * sin(X_u[3]);
    double zc = exp(X_u[1]) * cos(X_u[2]);

    double xs = sqrt(r_spot * r_spot + a * a) * cos(xspot[3]);
    double ys = sqrt(r_spot * r_spot + a * a) * sin(xspot[3]);
    double zs = r_spot * cos(xspot[2]);

    // distance^2 between photon position and spot center
    double xx = fabs(pow(xc - xs, 2) + pow(yc - ys, 2) + pow(zc - zs, 2));

    if (xx <= 4.)
        return exp(-(xx) / 2. / Rspot / Rspot);
    return 0.;
}

// Return emissivity for the thin disk model described in Dexter 2009
double emissivity_thindisk(double *X_u) {
    double r = logscale ? exp(X_u[1]) : X_u[1];
    return 1 / r / r;
}

// Planck function
double planck_function(double nu, double THETA_e) {
    double T = THETA_e * ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT /
               BOLTZMANN_CONSTANT;
    return 2. * PLANCK_CONSTANT * nu * nu * nu /
           (SPEED_OF_LIGHT * SPEED_OF_LIGHT) * 1. /
           (exp(PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T)) - 1.);
}

// Return absorption coefficient - assume local thermodynamical equilibrium so
// that Kirchoff's Law applies
double absorption_coeff_TH(double j_nu, double nu, double THETA_e) {
    double B_nu = planck_function(nu, THETA_e); // Planck function
    return j_nu / B_nu;
}
