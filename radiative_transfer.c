/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
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

// Return emission coefficient for thermal free-free radiation
double emission_coeff_FFTHERMAL(double nu, double n_e, double T) {
    double n_i = n_e; // Assume neutral hydrogen plasma
    double Z = 1.;    // For H, Z = 1
    double g_ff = 1.; // Gaunt factor

    double j_nu = 5.44e-39 * (Z * Z / sqrt(T)) * n_i * n_e * g_ff *
                  exp(-PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T));

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

double planck_function2(double nu, double Thetae) {
    double X = PLANCK_CONSTANT * nu /
               (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT * Thetae);

    if (X < 2.e-3)
        return ((2. * PLANCK_CONSTANT / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)) /
                (X / 24. * (24. + X * (12. + X * (4. + X)))));

    return ((2. * PLANCK_CONSTANT / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)) /
            (exp(X) - 1.));
}

// Return absorption coefficient - assume local thermodynamical equilibrium so
// that Kirchoff's Law applies
double absorption_coeff_TH(double j_nu, double nu, double THETA_e) {
    double B_nu = planck_function(nu, THETA_e); // Planck function
    return j_nu / B_nu;
}

// Compute the photon frequency in the plasma frame:
double freq_in_plasma_frame(double Uplasma_u[4], double k_d[4]) {
    double nu_plasmaframe = 0.;
    int i;
    LOOP_i nu_plasmaframe += Uplasma_u[i] * k_d[i];
    nu_plasmaframe *=
        -(ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / PLANCK_CONSTANT;

    if (isnan(nu_plasmaframe))
        fprintf(stderr, "NAN in plasma frame %e %e %e %e %e\n", nu_plasmaframe,
                Uplasma_u[0], Uplasma_u[1], Uplasma_u[2], Uplasma_u[3]);
    return nu_plasmaframe;
}

// See eqn 73 in Dexter 2016
double pitch_angle(double *X_u, double *k_u, double *B_u, double *Uplasma_u) {
    double b_dot_k = inner_product(X_u, B_u, k_u);
    double b_dot_b = inner_product(X_u, B_u, B_u);
    double k_dot_u = inner_product(X_u, k_u, Uplasma_u);

    // Compute and clamp result (result can slightly exceed domain of acos due
    // to numerics)
    double result = acos(b_dot_k / (-k_dot_u * sqrt(fabs(b_dot_b) + 1.e-15)));
    result = fmax(fmin(result, M_PI), 0.);

    //    return result;

    // NEW VERSION
    double B, k, mu;

    B = sqrt(fabs(inner_product(X_u, B_u, B_u)));

    if (B == 0.)
        return (M_PI / 2.);

    k = fabs(inner_product(X_u, k_u, Uplasma_u));

    mu = inner_product(X_u, k_u, B_u) / (k * B);

    if (fabs(mu) > 1.)
        mu /= fabs(mu);

    if (isnan(mu))
        fprintf(stderr, "isnan get_bk_angle\n");

    return (acos(mu));
}

/*
double radiative_transfer(double *lightpath, int steps, double frequency){
    int IN_VOLUME, path_counter;
    double I_current = 0.;
    double dI        = 0.;
    double j_nu      = 0.;
    double B, THETA_e, pitch_ang, nu_p, n_e, nu_p2, dl_current;
    int i;
    double X_u[4], k_u[4], k_d[4], B_u[4], Uplasma_u[4];
    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    double tau = 0.;

    double a_nu = 0.;

    double K_inv_old = 0, j_inv_old=0, dtau_old=0;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--){
        // Current position, wave vector, and dlambda
        LOOP_i{
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter-1) * 9 + 8]);

        // Obtain the parameters n_e, THETA_e, B, and Uplasma_u at X_u
        //get_plasma_parameters(X_u, &n_e, &THETA_e, &B, Uplasma_u);
        get_fluid_params(X_u, &n_e, &THETA_e, &B, B_u, Uplasma_u, &IN_VOLUME);

        // Check whether the ray is currently in the GRMHD simulation volume
        if(IN_VOLUME){
            // Obtain pitch angle: still no units (geometric)
            pitch_ang = pitch_angle(X_u, k_u, B_u, Uplasma_u);

            // CGS UNITS USED FROM HERE ON OUT
            //////////////////////////////////

            // Scale the wave vector to correct energy
            LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                             (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // Convert distance dlambda accordingly
            dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
(PLANCK_CONSTANT * frequency);

            //lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
            nu_p2 = nu_p * nu_p;

            // Obtain emission coefficient in current plasma conditions
           j_nu = emission_coeff_THSYNCHAV(B, THETA_e, nu_p, n_e);

            // Obtain absorption coefficient
            if (ABSORPTION){
                a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);
            }

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT / (ELECTRON_MASS * SPEED_OF_LIGHT *
SPEED_OF_LIGHT);

            double redshift = frequency / nu_p;

            double dtau  = (nu_p * a_nu * dl_current * C + dtau_old);
            double K_inv = (nu_p * a_nu);
            double j_inv = (j_nu / nu_p2);

            // Only add I_current if it is not NaN
            if(j_nu == j_nu && exp(X_u[1]) < RT_OUTER_CUTOFF){ // I_current +=
exp(-tau) * j_nu / nu_p / nu_p * dl_current * C;
       //         I_current += dI; // Old way of integrating
                double Ii=I_current;
                double S = j_inv/K_inv;
                if(K_inv == 0 )
                        I_current = Ii;
                else if(dtau < 1.e-5)
                        I_current = Ii - (Ii - S) * ( 0.166666667*dtau * (6. -
dtau * (3. - dtau))); else{ double efac = exp(-dtau); I_current = Ii*efac +
S*(1. - efac);
                }
             }

        }
    }

    // Store integrated intensity in the image
    return I_current * pow(frequency, 3.);
}
*/
