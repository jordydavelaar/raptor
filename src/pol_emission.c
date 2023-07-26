/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "gsl/gsl_sf_hyperg.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <gsl/gsl_sf_bessel.h>

// FUNCTIONS
////////////

// B.15:
double f_m(double X) {
    return 2.011 * exp(-pow(X, 1.035) / 4.7) -
           cos(X * 0.5) * exp(-pow(X, 1.2) / 2.73) - 0.011 * exp(-X / 47.2) +
           (0.011 * exp(-X / 47.2) - pow(2., -1. / 3.) / pow(3., 23. / 6.) *
                                         10000. * M_PI * pow(X, -8. / 3.)) *
               0.5 * (1. + tanh(10. * log(X / 120.)));
}

// Bessel function approximations:
double bessel_appr(int n, double x) {

    // use taylor expanded version
    if (x < 1. / 5.) {
        if (n == 0)
            return -log(x / 2.) - 0.5772;

        if (n == 1)
            return 1. / x;

        if (n == 2)
            return 2. / x / x;
    }
    // in this case all bessel functions are really small... Theta_e is small,
    // so no emission anyway.?
    else if (x > (1 / 0.00004) && 0) {
        return 1e-100;
    } else // use full versions in between.
        return gsl_sf_bessel_Kn(n, x);

    exit(0);
}

// Planck function
double planck_function(double nu, double THETA_e) {
    double T = THETA_e * ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT /
               BOLTZMANN_CONSTANT;
    return 2. * PLANCK_CONSTANT * nu * nu * nu /
           (SPEED_OF_LIGHT * SPEED_OF_LIGHT) * 1. /
           (exp(PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T)) - 1.);
}

double DeltaJ_5(double X) {
    return 0.43793091 * log(1. + 0.00185777 * pow(X, 1.50316886));
}
double get_w(double theta_e) {
    return (kappa-3)*theta_e/kappa;
    //return theta_e;
}
////////////
///////////////Rho_Q

double rho_Q_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B) {
    double w = get_w(theta_e);
    double factor_F = 0, factor_Q = 0;

    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));

    double F_factor_K3_5 =
        (1 - exp(-pow(X_kappa, 0.84) / 30.) -
         sin(X_kappa / 10.) * exp(-3 * pow(X_kappa, 0.471) / 2.));
    double F_factor_K4 = (1 - exp(-pow(X_kappa, 0.84) / 18.) -
                          sin(X_kappa / 6.) * exp(-7 * pow(X_kappa, 0.5) / 4.));
    double F_factor_K4_5 = (1 - exp(-pow(X_kappa, 0.84) / 12.) -
                            sin(X_kappa / 4.) * exp(-2 * pow(X_kappa, 0.525)));
    double F_factor_K5 =
        (1 - exp(-pow(X_kappa, 0.84) / 8.) -
         sin(3 * X_kappa / 8.) * exp(-9 * pow(X_kappa, 0.541) / 4.));

    double Q_factor_K3_5 =
        (17. * w - 3 * pow(w, 1. / 2.) + 7 * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K4 = ((46. / 3.) * w - (5. / 3.) * pow(w, 1. / 2.) +
                          (17. / 3.) * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K4_5 = (14. * w - (13. / 8.) * pow(w, 1. / 2.) +
                            (9. / 2.) * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K5 =
        ((25. / 2.) * w - pow(w, 1. / 2.) + 5 * pow(w, 1. / 2.) * exp(-5. * w));

    if (kappa == 3.5) {
        factor_F = F_factor_K3_5;
        ;
        factor_Q = Q_factor_K3_5;
    } else if (kappa == 4.0) {
        factor_F = F_factor_K4;
        factor_Q = Q_factor_K4;
    } else if (kappa == 4.5) {
        factor_F = F_factor_K4_5;
        factor_Q = Q_factor_K4_5;
    } else if (kappa == 5.0) {
        factor_F = F_factor_K5;
        factor_Q = Q_factor_K5;
    } else {
        fprintf(stderr, "kappa value not supported\n");
        exit(1);
    }
    return -n_e * pow(ELECTRON_CHARGE, 2) * pow(nuc, 2) * pow(sin(theta_B), 2) /
           (ELECTRON_MASS * SPEED_OF_LIGHT * pow(nu, 3)) * factor_F * factor_Q;
}

double rho_Q_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B) {
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    return 2. * M_PI * nu / 2. / SPEED_OF_LIGHT * wp2 * omega0 * omega0 /
           pow(2. * M_PI * nu, 4.) * f_m(Xe) *
           (bessel_appr(1, Thetaer) / bessel_appr(2, Thetaer) + 6. * theta_e) *
           sin(theta_B) * sin(theta_B);
}

double rho_Q_power(double theta_e, double n_e, double nu, double B,
                   double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double P_perp =
        (n_e * ELECTRON_CHARGE * ELECTRON_CHARGE) /
        (ELECTRON_MASS * SPEED_OF_LIGHT * nuc * sin(theta_B)) * (power - 1) *
        (1. / (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power)));
    return -P_perp * pow((nuc * sin(theta_B)) / nu, 3.) *
           (pow(gamma_min, 2. - power) / ((power / 2.) - 1.)) *
           (1. -
            pow((2. * nuc * sin(theta_B) * gamma_min * gamma_min) / (3. * nu),
                (power / 2.) - 1.));
}

double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
#if (DF == KAPPA)
    return rho_Q_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return rho_Q_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return rho_Q_thermal(theta_e, n_e, nu, B, theta_B);
#endif
}

///////////Rho_V

double rho_V_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B) {
    double w = get_w(theta_e);
    double factor_G = 0, factor_V = 0;
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2.) * sin(theta_B));
    double G_factor_K3_5 =
        (1. - 0.17 * log(1 + 0.447 * pow(X_kappa, -1. / 2.)));
    double G_factor_K4 = (1. - 0.17 * log(1 + 0.391 * pow(X_kappa, -1. / 2.)));
    double G_factor_K4_5 =
        (1. - 0.17 * log(1 + 0.348 * pow(X_kappa, -1. / 2.)));
    double G_factor_K5 = (1. - 0.17 * log(1 + 0.313 * pow(X_kappa, -1. / 2.)));

    double V_factor_K3_5 =
        ((pow(w, 2.) + 2. * w + 1.) / ((25. / 8.) * pow(w, 2.) + 4. * w + 1.));
    double V_factor_K4 = ((pow(w, 2.) + 54 * w + 50.) /
                          ((30. / 11.) * pow(w, 2) + 134 * w + 50));
    double V_factor_K4_5 = ((pow(w, 2.) + 43 * w + 38.) /
                            ((7. / 3.) * pow(w, 2) + (185. / 2.) * w + 38));
    double V_factor_K5 = ((w + (13. / 14.)) / (2 * w + (13. / 14.)));

    if (kappa == 3.5) {
        factor_G = G_factor_K3_5;
        ;
        factor_V = V_factor_K3_5;
    } else if (kappa == 4.0) {
        factor_G = G_factor_K4;
        factor_V = V_factor_K4;
    } else if (kappa == 4.5) {
        factor_G = G_factor_K4_5;
        factor_V = V_factor_K4_5;
    } else if (kappa == 5.0) {
        factor_G = G_factor_K5;
        factor_V = V_factor_K5;
    } else {
        fprintf(stderr, "kappa value not supported\n");
        exit(1);
    }

    return 2 * n_e * pow(ELECTRON_CHARGE, 2.) * nuc * cos(theta_B) /
           (ELECTRON_MASS * SPEED_OF_LIGHT * pow(nu, 2.)) *
           bessel_appr(0., 1. / w) / (bessel_appr(2., 1. / w)) * factor_G *
           factor_V;
}

double rho_V_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B) {
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    double fit_factor;
    double k2 = bessel_appr(2, Thetaer);
    double k0 = bessel_appr(0, Thetaer);

#if(DEXTER)
    fit_factor = (k0 - DeltaJ_5(Xe)) / k2;
#else
    double shgmfunc = 1 - 0.11*log(1 + 0.035*Xe);
    double k_ratio = (k2 > 0) ? k0/k2 : 1;

    fit_factor = k_ratio * shgmfunc;
#endif

    return 2.0 * M_PI * nu / SPEED_OF_LIGHT * wp2 * omega0 /
           pow(2. * M_PI * nu, 3.) * fit_factor * cos(theta_B);

}

double rho_V_power(double theta_e, double n_e, double nu, double B,
                   double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double P_perp =
        (n_e * ELECTRON_CHARGE * ELECTRON_CHARGE) /
        (ELECTRON_MASS * SPEED_OF_LIGHT * nuc * sin(theta_B)) * (power - 1.) *
        pow((pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power)), -1.);
    return 2. * P_perp * ((power + 2.) / (power + 1.)) *
           pow(((nuc * sin(theta_B)) / (nu)), 2.) *
           pow(gamma_min, -(power + 1.)) * log(gamma_min) * cos(theta_B) /
           sin(theta_B);
}
double rho_V(double theta_e, double n_e, double nu, double B, double theta_B) {
#if (DF == KAPPA)
    return rho_V_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return rho_V_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return rho_V_thermal(theta_e, n_e, nu, B, theta_B);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////J
/// for thermal and kappa
double j_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_I = 1;
    double J_high_factor_I = 1;
    double J_S_factor_I = 1;
    double J_x_I = 3 * pow(kappa, -3. / 2.);
    double J_I_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_I;
    double J_I_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_I;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac *
           pow((pow(J_I_low_kappa, -J_x_I) + pow(J_I_high_kappa, -J_x_I)),
               -1. / J_x_I) *
           J_S_factor_I;
}

double I_I(double x) {
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) + 0.9977 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double j_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                  (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                  theta_e;

    double x = nu / nu_c;

    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_I(x);
}

double j_I_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_I_power_factor = 1;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac * (pow(3, power / 2.) * (power - 1) * sin(theta_B)) /
           (2 * (power + 1) *
            (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power - 1) / 12.) * tgamma((3 * power + 19) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power - 1) / 2.) *
           j_I_power_factor;
}

double j_I(double theta_e, double n_e, double nu, double B, double theta_B) {
#if (DF == KAPPA)
    return j_I_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return j_I_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return j_I_thermal(theta_e, n_e, nu, B, theta_B);
#endif
}
double j_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_Q = 0.5;
    double J_high_factor_Q = (pow((4. / 5.), 2.) + kappa / 50.);
    double J_S_factor_Q = (-1);
    double J_x_Q = (37. / 10.) * pow(kappa, -8. / 5.);
    double J_Q_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_Q;
    double J_Q_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_Q;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    //this needs to  change sign in our IEEE convention! 
    return - prefac *
           pow((pow(J_Q_low_kappa, -J_x_Q) + pow(J_Q_high_kappa, -J_x_Q)),
               -1. / J_x_Q) *
           J_S_factor_Q;
}

double I_Q(double x) {
    return 2.5651 *
           (1 + 0.93193 * pow(x, -1. / 3.) + 0.499873 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double j_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                  (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                  theta_e;
    double x = nu / nu_c;
    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_Q(x);
}
double j_Q_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_Q_power_factor = -(power + 1) / (power + 7. / 3.);
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;

    //this needs to  change sign in our IEEE convention! 
    return - prefac * (pow(3, power / 2.) * (power - 1) * sin(theta_B)) /
           (2 * (power + 1) *
            (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power - 1) / 12.) * tgamma((3 * power + 19) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power - 1) / 2.) *
           j_Q_power_factor;
}

double j_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
#if (DF == KAPPA)
    return j_Q_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return j_Q_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return j_Q_thermal(theta_e, n_e, nu, B, theta_B);
#endif
}

double j_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_V =
        (pow(0.75, 2.) * pow((pow(sin(theta_B), -12. / 5.) - 1), 12. / 25.) *
         pow(kappa, -66. / 125.) * pow(X_kappa, -7. / 20.) / w);
    double J_high_factor_V =
        (pow((7. / 8.), 2.) *
         pow((pow(sin(theta_B), -5. / 2.) - 1), 11. / 25.) *
         pow(kappa, -11. / 25.) * pow(X_kappa, -1. / 2.) / w);
    double J_S_factor_V = (sign(cos(theta_B)));
    double J_x_V = 3 * pow(kappa, -3. / 2.);
    double J_V_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_V;
    double J_V_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_V;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac *
           pow((pow(J_V_low_kappa, -J_x_V) + pow(J_V_high_kappa, -J_x_V)),
               -1. / J_x_V) *
           J_S_factor_V;
}

double I_V(double x) {
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
            0.0292545 * pow(x, -0.5) + 2.03773 * pow(x, -1. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double j_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                  (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                  theta_e;
    double x = nu / nu_c;
    return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / tan(theta_B) /
           3. / sqrt(3.) / SPEED_OF_LIGHT / theta_e / theta_e / theta_e *
           I_V(x);
}

double j_V_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_V_power_factor = (171. / 250.) *
                              (pow(power, 49. / 100.) / tan(theta_B)) *
                              pow(nu / (3. * nuc * sin(theta_B)), -1. / 2.);
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac * (pow(3., power / 2.) * (power - 1.) * sin(theta_B)) /
           (2. * (power + 1.) *
            (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power))) *
           tgamma((3. * power - 1.) / 12.) * tgamma((3. * power + 19.) / 12.) *
           pow(nu / (nuc * sin(theta_B)), -(power - 1.) / 2.) *
           j_V_power_factor;
}

double j_V(double theta_e, double n_e, double nu, double B, double theta_B) {
#if (DF == KAPPA)
    return j_V_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return j_V_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return j_V_thermal(theta_e, n_e, nu, B, theta_B);
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////A
/// for thermal and kappa
double hyp2F1_f(double theta_e) {
    double w = get_w(theta_e);
    double a = kappa - 1. / 3.;
    double b = kappa + 1.;
    double c = kappa + 2. / 3.;
    double X = kappa * w;
    double z;

    if (X < 1e-4) {
        return 0;
    }

    else {
        z = -X;
        return pow(1. - z, -a) * tgamma(c) * tgamma(b - a) /
                   (tgamma(b) * tgamma(c - a)) *
                   gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z)) +
               pow(1. - z, -b) * tgamma(c) * tgamma(a - b) /
                   (tgamma(a) * tgamma(c - b)) *
                   gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1. - z));
    }

    return 0;
}

double a_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_I = 1;
    double A_high_factor_I = (pow((3.0 / kappa), 19. / 4.) + (3. / 5.));
    double A_S_factor_I = 1;
    double A_x_I = pow((-7. / 4. + kappa * 8. / 5.), -43. / 50.);
    double A_I_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e) * A_low_factor_I;
    double A_I_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_I;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac *
           pow((pow(A_I_low_kappa, -A_x_I) + pow(A_I_high_kappa, -A_x_I)),
               -1. / A_x_I) *
           A_S_factor_I;
}

double a_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_I_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_I_thermal / B_nu;
}

double a_I_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_I_power_factor = 1;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac * (pow(3, (power + 1) / 2.) * (power - 1)) /
           (4 * (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power + 2) / 12.) * tgamma((3 * power + 22) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2) / 2.) *
           a_I_power_factor;
}

double a_I(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_I_thermal) {
#if (DF == KAPPA)
    return a_I_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return a_I_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return a_I_thermal(theta_e, n_e, nu, B, theta_B, j_I_thermal);
#endif
}

double a_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_Q = (25. / 48.);
    double A_high_factor_Q =
        (pow(21, 2.) * pow(kappa, -pow(12. / 5., 2.)) + (11. / 20.));
    double A_S_factor_Q = (-1);
    double A_x_Q = (7. / 5.) * pow(kappa, -23. / 20.);
    double A_Q_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e) * A_low_factor_Q;
    double A_Q_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_Q;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);

    //this needs to  change sign in our IEEE convention! 
    return -prefac *
           pow((pow(A_Q_low_kappa, -A_x_Q) + pow(A_Q_high_kappa, -A_x_Q)),
               -1. / A_x_Q) *
           A_S_factor_Q;
}

double a_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_Q_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_Q_thermal / B_nu;
}

double a_Q_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_Q_power_factor =
        -pow(((17 * power / 500.) - 43. / 1250.), 43. / 500.);
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);

    //this needs to  change sign in our IEEE convention! 
    return -prefac * (pow(3, (power + 1) / 2.) * (power - 1)) /
           (4 * (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power + 2) / 12.) * tgamma((3 * power + 22) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2) / 2.) *
           a_Q_power_factor;
}

double a_Q(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_Q_thermal) {
#if (DF == KAPPA)
    return a_Q_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return a_Q_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return a_Q_thermal(theta_e, n_e, nu, B, theta_B, j_Q_thermal);
#endif
}

double a_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double w = get_w(theta_e);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_V =
        ((77. / (100 * w)) *
         pow((pow(sin(theta_B), -114. / 50.) - 1), 223. / 500.) *
         pow(X_kappa, -7. / 20.) * pow(kappa, -7. / 10.));
    double A_high_factor_V =
        ((143. / 10.) * pow(w, -116. / 125.) *
         pow((pow(sin(theta_B), -41. / 20.) - 1), 1. / 2.) *
         (pow(13, 2) * pow(kappa, -8) + (13. / 2500.) * kappa - (263. / 5000.) +
          (47. / (200 * kappa))) *
         pow(X_kappa, -1. / 2.));
    double A_S_factor_V = (sign(cos(theta_B)));
    double A_x_V = ((61. / 50.) * pow(kappa, -142. / 125.) + (7. / 1000.));
    double A_V_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e) * A_low_factor_V;
    double A_V_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_V;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);

//from ipole/symphony. 
        /*The Stokes V absorption coefficient changes sign at observer_angle
	    equals 90deg, but this formula does not.  This discrepancy is a
	    bug in this formula, and is patched by the term below.*/
	  double sign_bug_patch = cos(theta_B)/fabs(cos(theta_B));

	  /*NOTE: Sign corrected; the sign in Leung et al. (2011)
	    and Pandya et al. (2016) for Stokes V transfer coefficients
	    does not follow the convention the papers describe (IEEE/IAU);
	    the sign has been corrected here.*/

    return -sign_bug_patch * prefac *
           pow((pow(A_V_low_kappa, -A_x_V) + pow(A_V_high_kappa, -A_x_V)),
               -1. / A_x_V) *
           A_S_factor_V;
}

double a_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_V_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_V_thermal / B_nu;
}

double a_V_power(double theta_e, double n_e, double nu, double B,
                 double theta_B) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_V_power_factor =
        pow(((71. * power / 100.) + 22. / 625.), 197. / 500.) *
        pow(((31. / 10.) * pow(sin(theta_B), -48. / 25.) - 31. / 10.),
            64. / 125.) *
        pow(nu / (nuc * sin(theta_B)), -1. / 2.) * sign(cos(theta_B));
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac * (pow(3., (power + 1) / 2.) * (power - 1.)) /
           (4. * (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power))) *
           tgamma((3. * power + 2.) / 12.) * tgamma((3. * power + 22.) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2.) / 2.) *
           a_V_power_factor;
}

double a_V(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_V_thermal) {
#if (DF == KAPPA)
    return a_V_kappa(theta_e, n_e, nu, B, theta_B);
#elif (DF == POWER)
    return a_V_power(theta_e, n_e, nu, B, theta_B);
#elif (DF == TH)
    return a_V_thermal(theta_e, n_e, nu, B, theta_B, j_V_thermal);
#endif
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
