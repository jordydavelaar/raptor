/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "functions.h"
#include "parameters.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_parallel(const double y[], double complex f_u[], double fvector[],
                double complex f_u_vector[]) {
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];
    // Einstein summation over indices v and w

    LOOP_ijk gamma_udd[i][j][k] = 0.;

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double complex A_u[4] = {0., 0., 0., 0.}; // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
#if (metriic == MKSBHAC || metric == MKSHARM)
    connection_udd(X_u, gamma_udd);
#else
    connection_num_udd(X_u, gamma_udd);
#endif

    // Compute 4-acceleration using the geodesic equation
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];
    LOOP_i {
        fvector[i] = U_u[i];
        fvector[i + 4] = A_u[i];
    }

    // Reset A_u
    LOOP_i A_u[i] = 0.;

    // Compute f_u vector acceleration
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * f_u[k];
    LOOP_i { f_u_vector[i] = A_u[i]; }
}

void rk4_step_f(double y[], double complex f_u[], double dt) {
    // Array containing all "update elements" (4 times Nelements because RK4)
    double dx[4 * 2 * 4];
    double complex df[4 * 4];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK4
    double yshift[4 * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};
    double complex f_u_shift[4] = {f_u[0], f_u[1], f_u[2], f_u[3]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[4 * 2];
    double complex f_u_vector[4];

    // Compute the RK4 update coefficients ('K_n' in lit., 'dx' here)
    int i, q;
    double complex weights[4] = {0.5, 0.5, 1.,
                                 0.}; // Weights used for updating y
    for (q = 0; q < 4; q++) {
        f_parallel(
            yshift, f_u_shift, fvector,
            f_u_vector); // Apply function f to current y to obtain fvector
        for (i = 0; i < 4 * 2; i++) {
            dx[q * 4 * 2 + i] = dt * fvector[i]; // Use fvector to fill dx
            yshift[i] = y[i] + dx[q * 4 * 2 + i] * weights[q]; // Update y
        }
        for (i = 0; i < 4; i++) {
            df[q * 4 + i] = dt * f_u_vector[i];
            f_u_shift[i] = f_u[i] + df[q * 4 + i] * weights[q];
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < 4 * 2; i++) {
        y[i] = y[i] + 1. / 6. *
                          (dx[0 * 4 * 2 + i] + dx[1 * 4 * 2 + i] * 2. +
                           dx[2 * 4 * 2 + i] * 2. + dx[3 * 4 * 2 + i]);
    }

    // Update the f-vector (polarization)
    for (i = 0; i < 4; i++) {
        f_u[i] = f_u[i] + 1. / 6. *
                              (df[0 * 4 + i] + df[1 * 4 + i] * 2. +
                               df[2 * 4 + i] * 2. + df[3 * 4 + i]);
    }
}

void f_tetrad_to_stokes(double Iinv, double Iinv_pol,
                        double complex f_tetrad_u[], double complex S_A[]) {
    S_A[0] = Iinv;
    S_A[1] = Iinv_pol * (cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) -
                         cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]));
    S_A[2] = Iinv_pol * (conj(f_tetrad_u[1]) * f_tetrad_u[2] +
                         f_tetrad_u[1] * conj(f_tetrad_u[2]));
    S_A[3] = Iinv_pol * (I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] -
                              f_tetrad_u[1] * conj(f_tetrad_u[2])));
}

void stokes_to_f_tetrad(double complex S_A[], double *Iinv, double *Iinv_pol,
                        double complex f_tetrad_u[]) {

    *Iinv = S_A[0];

    *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

    double Qnorm = S_A[1] / (*Iinv_pol);
    double Unorm = S_A[2] / (*Iinv_pol);
    double Vnorm = S_A[3] / (*Iinv_pol);

    // source:
    // https://physics.stackexchange.com/questions/238957/converting-stokes-parameters-to-jones-vector
    f_tetrad_u[1] = sqrt((1. + Qnorm) / 2.);

    if (f_tetrad_u[1] == 0)
        f_tetrad_u[2] = 1.;
    else
        f_tetrad_u[2] =
            Unorm / (2. * f_tetrad_u[1]) - I * Vnorm / (2. * f_tetrad_u[1]);
}

// NOTE: works only in Kerr metric
// Ziri's suggestion: construct U vecs
void construct_U_vector(const double X_u[], double U_u[]) {
// Obtain relevant metric terms:
#if (metric == CKS)
    double U_KS[4];
    double X_KS[4];

    CKS_to_KS(X_u, X_KS);

    double g_uu[4][4];
    metric_KS_uu(X_KS, g_uu);
    double g_uu00 = g_uu[0][0];
    double g_uu03 = g_uu[0][3];
    double g_uu33 = g_uu[3][3];

#else
    double g_uu[4][4];
    metric_uu(X_u, g_uu);
    double g_uu00 = g_uu[0][0];
    double g_uu03 = g_uu[0][3];
    double g_uu33 = g_uu[3][3];
#endif
    // Observer/plasma wave vector:
    double U_d[4] = {-1., 0., 0., 0.};
    double B__ = -g_uu03 * U_d[0] / g_uu33;
    double C__ = -(1. + g_uu00 * U_d[0] * U_d[0]) / g_uu33;

    // Properly normalize U_u:
    U_d[3] = B__ + sqrt(B__ * B__ + C__);

#if (metric == CKS)
    LOOP_i {
        U_KS[i] = 0.;
        U_u[i] = 0;
    }
    raise_index_KS(X_KS, U_d, U_KS);

    double coordKS[8];
    double coordCKS[8];

    LOOP_i {
        coordKS[i] = X_KS[i];
        coordKS[i + 4] = U_KS[i];
    }

    KS_to_CKS_u(coordKS, coordCKS);

    LOOP_i U_u[i] = coordCKS[i + 4];

#else
    LOOP_i U_u[i] = 0.;
    raise_index(X_u, U_d, U_u);

#endif
}

// NEW FUNCTIONS JUNE 2021
//////////////////////////

// Transform f_tetrad_u to f_u
void f_tetrad_to_f(double complex *f_u, double tetrad_u[][4],
                   double complex *f_tetrad_u) {

    LOOP_i f_u[i] = 0.;
    LOOP_ij f_u[i] += tetrad_u[i][j] * f_tetrad_u[j];
}

// Transform f_u to f_tetrad_u
void f_to_f_tetrad(double complex *f_tetrad_u, double tetrad_d[][4],
                   double complex *f_u) {

    LOOP_i f_tetrad_u[i] = 0.;
    LOOP_ij f_tetrad_u[i] += tetrad_d[j][i] * f_u[j];
}

void evaluate_coeffs_user(double *jI, double *jQ, double *jU, double *jV,
                          double *rQ, double *rU, double *rV, double *aI,
                          double *aQ, double *aU, double *aV, double nu_p,
                          struct GRMHD modvar, double pitch_ang) {
    double jI_thermal, jQ_thermal, jV_thermal, jU_thermal = 0;
    double jI_kappa, jQ_kappa, jV_kappa, jU_kappa = 0;

    double aI_kappa, aV_kappa, aQ_kappa, aU_kappa = 0;
    double aI_thermal, aV_thermal, aQ_thermal, aU_thermal = 0;

    double eps, epsilon;

    jI_thermal =
        j_I_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    jV_thermal =
        j_V_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    jQ_thermal =
        j_Q_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

    aI_thermal = a_I_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B,
                             pitch_ang, jI_thermal);
    aV_thermal = a_V_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B,
                             pitch_ang, jV_thermal);
    aQ_thermal = a_Q_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B,
                             pitch_ang, jQ_thermal);

    jI_kappa = j_I_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    jV_kappa = j_V_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    jQ_kappa = j_Q_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

    aI_kappa = a_I_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    aV_kappa = a_V_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    aQ_kappa = a_Q_kappa(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

    // to invariant forms...
    jI_thermal /= (nu_p * nu_p);
    jV_thermal /= (nu_p * nu_p);
    jQ_thermal /= (nu_p * nu_p);

    aI_thermal *= nu_p;
    aV_thermal *= nu_p;
    aQ_thermal *= nu_p;

    jI_kappa /= (nu_p * nu_p);
    jV_kappa /= (nu_p * nu_p);
    jQ_kappa /= (nu_p * nu_p);

    aI_kappa *= nu_p;
    aV_kappa *= nu_p;
    aQ_kappa *= nu_p;

    // MIXED MODEL
    epsilon = 1.;
    eps = epsilon * (1. - exp(-pow(modvar.beta, -2.))) *
          (1 - exp(-pow(modvar.sigma / modvar.sigma_min, 2)));

    *jI = (1. - eps) * jI_thermal + eps * jI_kappa;
    *jV = (1. - eps) * jV_thermal + eps * jV_kappa;
    *jU = 0.0;
    *jQ = (1. - eps) * jQ_thermal + eps * jQ_kappa;

    *aV = (1. - eps) * aV_thermal + eps * aV_kappa;
    *aQ = (1. - eps) * aQ_thermal + eps * aQ_kappa;
    *aU = 0.0;
    *aI = (1. - eps) * aI_thermal + eps * aI_kappa;
}

void evaluate_coeffs_single(double *jI, double *jQ, double *jU, double *jV,
                            double *rQ, double *rU, double *rV, double *aI,
                            double *aQ, double *aU, double *aV, double nu_p,
                            struct GRMHD modvar, double pitch_ang) {
    *jI = j_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    *jQ = j_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    *jU = 0.;
    *jV = j_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

    *rQ = rho_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
    *rU = 0.;
    *rV = rho_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

    *aI = a_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jI);
    *aQ = a_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jQ);
    *aU = 0;
    *aV = a_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jV);

    // Transform to invariant forms
    *jI /= (nu_p * nu_p);
    *jQ /= (nu_p * nu_p);
    *jV /= (nu_p * nu_p);

    *aI *= nu_p;
    *aQ *= nu_p;
    *aV *= nu_p;

    *rQ *= nu_p;
    *rV *= nu_p;
}
int check_stiffness(double jI, double jQ, double jU, double jV, double rQ,
                    double rU, double rV, double aI, double aQ, double aU,
                    double aV, double dl_current) {
    // int STIFF = check_stiffness...
    double a2 = rQ * rQ + rV * rV - aQ * aQ - aV * aV;
    double a0 = -2. * aV * aQ * rV * rQ - aQ * aQ * rQ * rQ - aV * aV * rV * rV;

    complex double zplus = (-a2 + sqrt(a2 * a2 - 4. * a0)) / 2.;
    complex double zminus = (-a2 - sqrt(a2 * a2 - 4. * a0)) / 2.;

    complex double l1 = aI + sqrt(zplus);
    complex double l2 = aI - sqrt(zplus);
    complex double l3 = aI + sqrt(zminus);
    complex double l4 = aI - sqrt(zminus);

    complex double tau1 = dl_current * l1;
    complex double tau2 = dl_current * l2;
    complex double tau3 = dl_current * l3;
    complex double tau4 = dl_current * l4;

    complex double mag1 = 1. + tau1 + 0.5 * tau1 * tau1 +
                          1. / 6. * tau1 * tau1 * tau1 +
                          1. / 24. * tau1 * tau1 * tau1 * tau1;
    complex double mag2 = 1. + tau2 + 0.5 * tau2 * tau2 +
                          1. / 6. * tau2 * tau2 * tau2 +
                          1. / 24. * tau2 * tau2 * tau2 * tau2;
    complex double mag3 = 1. + tau3 + 0.5 * tau3 * tau3 +
                          1. / 6. * tau3 * tau3 * tau3 +
                          1. / 24. * tau3 * tau3 * tau3 * tau3;
    complex double mag4 = 1. + tau4 + 0.5 * tau4 * tau4 +
                          1. / 6. * tau4 * tau4 * tau4 +
                          1. / 24. * tau4 * tau4 * tau4 * tau4;

    double res1 = sqrt(mag1 * conj(mag1));
    double res2 = sqrt(mag2 * conj(mag2));
    double res3 = sqrt(mag3 * conj(mag3));
    double res4 = sqrt(mag4 * conj(mag4));

    int STIFF = 0;

    double STIFFTHRESH = 0.99;

    if (res1 > STIFFTHRESH || res2 > STIFFTHRESH || res3 > STIFFTHRESH ||
        res4 > STIFFTHRESH)
        STIFF = 1;

    return STIFF;
}

void pol_rte_rk4_step(double jI, double jQ, double jU, double jV, double rQ,
                      double rU, double rV, double aI, double aQ, double aU,
                      double aV, double dl_current, double C,
                      double complex S_A[]) {
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    // RK4 with constant coefficients
    // k1
    double complex Ik1 =
        dl_current * C * jI -
        dl_current * C * (aI * I0 + aQ * Q0 + aU * U0 + aV * V0);
    double complex Qk1 =
        dl_current * C * jQ -
        dl_current * C * (aQ * I0 + aI * Q0 + rV * U0 - rU * V0);
    double complex Uk1 =
        dl_current * C * jU -
        dl_current * C * (aU * I0 - rV * Q0 + aI * U0 + rQ * V0);
    double complex Vk1 =
        dl_current * C * jV -
        dl_current * C * (aV * I0 + rU * Q0 - rQ * U0 + aI * V0);

    // k2
    double complex Ik2 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik1) + aQ * (Q0 + 0.5 * Qk1) +
                              aU * (U0 + 0.5 * Uk1) + aV * (V0 + 0.5 * Vk1));
    double complex Qk2 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik1) + aI * (Q0 + 0.5 * Qk1) +
                              rV * (U0 + 0.5 * Uk1) - rU * (V0 + 0.5 * Vk1));
    double complex Uk2 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik1) - rV * (Q0 + 0.5 * Qk1) +
                              aI * (U0 + 0.5 * Uk1) + rQ * (V0 + 0.5 * Vk1));
    double complex Vk2 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik1) + rU * (Q0 + 0.5 * Qk1) -
                              rQ * (U0 + 0.5 * Uk1) + aI * (V0 + 0.5 * Vk1));

    // k3
    double complex Ik3 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik2) + aQ * (Q0 + 0.5 * Qk2) +
                              aU * (U0 + 0.5 * Uk2) + aV * (V0 + 0.5 * Vk2));
    double complex Qk3 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik2) + aI * (Q0 + 0.5 * Qk2) +
                              rV * (U0 + 0.5 * Uk2) - rU * (V0 + 0.5 * Vk2));
    double complex Uk3 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik2) - rV * (Q0 + 0.5 * Qk2) +
                              aI * (U0 + 0.5 * Uk2) + rQ * (V0 + 0.5 * Vk2));
    double complex Vk3 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik2) + rU * (Q0 + 0.5 * Qk2) -
                              rQ * (U0 + 0.5 * Uk2) + aI * (V0 + 0.5 * Vk2));

    // k4
    double complex Ik4 =
        dl_current * C * jI - dl_current * C *
                                  (aI * (I0 + Ik3) + aQ * (Q0 + Qk3) +
                                   aU * (U0 + Uk3) + aV * (V0 + Vk3));
    double complex Qk4 =
        dl_current * C * jQ - dl_current * C *
                                  (aQ * (I0 + Ik3) + aI * (Q0 + Qk3) +
                                   rV * (U0 + Uk3) - rU * (V0 + Vk3));
    double complex Uk4 =
        dl_current * C * jU - dl_current * C *
                                  (aU * (I0 + Ik3) - rV * (Q0 + Qk3) +
                                   aI * (U0 + Uk3) + rQ * (V0 + Vk3));
    double complex Vk4 =
        dl_current * C * jV - dl_current * C *
                                  (aV * (I0 + Ik3) + rU * (Q0 + Qk3) -
                                   rQ * (U0 + Uk3) + aI * (V0 + Vk3));

    S_A[0] = I0 + 1. / 6. * (Ik1 + 2. * Ik2 + 2. * Ik3 + Ik4);
    S_A[1] = Q0 + 1. / 6. * (Qk1 + 2. * Qk2 + 2. * Qk3 + Qk4);
    S_A[2] = U0 + 1. / 6. * (Uk1 + 2. * Uk2 + 2. * Uk3 + Uk4);
    S_A[3] = V0 + 1. / 6. * (Vk1 + 2. * Vk2 + 2. * Vk3 + Vk4);
}

void pol_rte_trapezoid_step(double jI, double jQ, double jU, double jV,
                            double rQ, double rU, double rV, double aI,
                            double aQ, double aU, double aV, double dl_current,
                            double C, double complex S_A[]) {
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    double u11 = 1. + 0.5 * dl_current * C * aI;
    double u12 = 0.5 * dl_current * C * aQ;
    double u14 = 0.5 * dl_current * C * aV;
    double l21 = 0.5 * dl_current * C * aQ / u11;
    double u22 = 1. + 0.5 * dl_current * C * aI - l21 * u12;
    double u23 = 0.5 * dl_current * C * rV;
    double u24 = -l21 * u14;
    double l32 = -0.5 * dl_current * C * rV / u22;
    double u33 = 1. + 0.5 * dl_current * C * aI - l32 * u23;
    double u34 = 0.5 * dl_current * C * rQ - l32 * u24;
    double l41 = 0.5 * dl_current * C * aV / u11;
    double l42 = -l41 * u12 / u22;
    double l43 = (-0.5 * dl_current * C * rQ - l42 * u23) / u33;
    double u44 =
        1. + 0.5 * dl_current * C * aI - l41 * u14 - l42 * u24 - l43 * u34;

    // Construct b-vector.
    double b1 =
        I0 + dl_current * C / 2. * (2. * jI - (aI * I0 + aQ * Q0 + aV * V0));
    double b2 =
        Q0 + dl_current * C / 2. * (2. * jQ - (aQ * I0 + aI * Q0 + rV * U0));
    double b3 =
        U0 + dl_current * C / 2. * (2. * jU - (-rV * Q0 + aI * U0 + rQ * V0));
    double b4 =
        V0 + dl_current * C / 2. * (2. * jV - (aV * I0 - rQ * U0 + aI * V0));

    // Construct y.
    double y1 = b1;
    double y2 = b2 - l21 * y1;
    double y3 = b3 - l32 * y2;
    double y4 = b4 - l41 * y1 - l42 * y2 - l43 * y3;

    // Construct x.
    double x4 = y4 / u44;
    double x3 = (y3 - u34 * x4) / u33;
    double x2 = (y2 - u23 * x3 - u24 * x4) / u22;
    double x1 = (y1 - u12 * x2 - u14 * x4) / u11;

    S_A[0] = x1;
    S_A[1] = x2;
    S_A[2] = x3;
    S_A[3] = x4;
}

void f_to_stokes(double complex f_u[], double complex f_tetrad_u[],
                 double tetrad_d[][4], double complex S_A[], double Iinv,
                 double Iinv_pol) {
    f_to_f_tetrad(f_tetrad_u, tetrad_d, f_u);

    // Get Stokes params from f_tetrad_u
    f_tetrad_to_stokes(Iinv, Iinv_pol, f_tetrad_u, S_A);
}

void stokes_to_f(double complex f_u[], double complex f_tetrad_u[],
                 double tetrad_u[][4], double complex S_A[], double *Iinv,
                 double *Iinv_pol) {
    stokes_to_f_tetrad(S_A, Iinv, Iinv_pol, f_tetrad_u);

    f_tetrad_to_f(f_u, tetrad_u, f_tetrad_u);
}

void pol_integration_step(struct GRMHD modvar, double frequency,
                          double *dl_current, double C, double X_u[],
                          double k_u[], double k_d[], int *POLARIZATION_ACTIVE,
                          double complex f_u[], double complex f_tetrad_u[],
                          double tetrad_d[][4], double tetrad_u[][4],
                          double complex S_A[], double *Iinv,
                          double *Iinv_pol) {

    double jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV;
    double pitch_ang, nu_p;
    // Unpolarized: 1) Create light path by integration. 2) For each
    // step in lightpath, perform one radiative transfer step.
    // Polarized:   1) Create light path by integration. 2) For each
    // step in lightpath, perform one radiative transfer step, AND,
    // OUTSIDE in_volume loop, do a spacetime propagation step.

    // TRANSFER STEP
    ////////////////

    // Obtain pitch angle: still no units (geometric)
    pitch_ang = pitch_angle(X_u, k_u, modvar.B_u, modvar.U_u);

    // CGS UNITS USED FROM HERE ON OUT
    //////////////////////////////////

    // Scale the wave vector to correct energy
    LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                     (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    // Convert distance dlambda accordingly
    *dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                   (PLANCK_CONSTANT * frequency);

    // lower the index of the wavevector
    lower_index(X_u, k_u, k_d);

    // Compute the photon frequency in the plasma frame:
    nu_p = freq_in_plasma_frame(modvar.U_u, k_d);

    // POLARIZED EMISSION/ABSORPTION COEFFS
    ///////////////////////////////////////
    evaluate_coeffs_single(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI, &aQ, &aU,
                           &aV, nu_p, modvar, pitch_ang);

    // Create tetrad, needed whether POLARIZATION_ACTIVE is true or
    // false.
    create_observer_tetrad(X_u, k_u, modvar.U_u, modvar.B_u, tetrad_u);
    create_tetrad_d(X_u, tetrad_u, tetrad_d);

    // FROM F VECTOR TO STOKES (when applicable)
    ////////////////////////////////////////////

    // If (POLARIZATION_ACTIVE), get Stokes params from f_u and p.
    // (Otherwise, never been in volume before; we simply use
    // S_I_current)
    if (*POLARIZATION_ACTIVE) {
        f_to_stokes(f_u, f_tetrad_u, tetrad_d, S_A, *Iinv, *Iinv_pol);
    }

    // Given Stokes params and plasma coeffs, compute NEW Stokes params
    // after plasma step.

    int STIFF = check_stiffness(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                                *dl_current);

    // If both rotation coeffs (times dlambda) are smaller than
    // threshold, take an RK4 step; otherwise, implicit Euler.
    // if (fabs(rQ) < THRESH && fabs(rV) < THRESH) {
    if (!STIFF) {
        pol_rte_rk4_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                         *dl_current, C, S_A);
    } else {
        pol_rte_trapezoid_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                               *dl_current, C, S_A);
    }

    // FROM STOKES TO F VECTOR
    ///////////////////////////

    *Iinv = S_A[0];
    *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

    //        fprintf(stderr,"Iinv %e Iinv_pol %e\n",*Iinv,*Iinv_pol);
    //        fprintf(stderr,"jI %e jQ %e jU %e jV %e\n",jI,jQ,jU,jV);
    //        check_tetrad_identities(X_u, tetrad_u);
    //        check_tetrad_compact(X_u, tetrad_u);

    // We have now updated the Stokes vector using plasma at current
    // position. Only do stuff below this line IF S_A[0] > 1.e-40. If
    // not, POLARIZATION_ACTIVE is set to FALSE and we reset S_A[i] = 0
    if (*Iinv_pol > 1.e-100) {
        stokes_to_f(f_u, f_tetrad_u, tetrad_u, S_A, Iinv, Iinv_pol);

        // Set POLARIZATION_ACTIVE to true; we are, after all,
        // in_volume.
        *POLARIZATION_ACTIVE = 1;

    } else {
        *POLARIZATION_ACTIVE = 0;
        S_A[1] = 0.;
        S_A[2] = 0.;
        S_A[3] = 0.;
    }
}

void construct_f_obs_tetrad_u(double *X_u, double *k_u, double complex *f_u,
                              double complex *f_obs_tetrad_u) {

    double cam_up_u[4] = {0., 0., 0., -1.};
    double U_obs_u[4] = {0., 0., 0., 0.};
    double obs_tetrad_u[4][4], obs_tetrad_d[4][4];
    LOOP_ij obs_tetrad_u[i][j] = 0.;
    LOOP_ij obs_tetrad_d[i][j] = 0.;

    construct_U_vector(X_u, U_obs_u);
    create_observer_tetrad(X_u, k_u, U_obs_u, cam_up_u, obs_tetrad_u);
    create_tetrad_d(X_u, obs_tetrad_u, obs_tetrad_d);

    // Convert f_u to f_obs_tetrad_u
    LOOP_i f_obs_tetrad_u[i] = 0.;
    LOOP_ij f_obs_tetrad_u[i] += obs_tetrad_d[j][i] * f_u[j];
}

void radiative_transfer_polarized(double *lightpath, int steps,
                                  double frequency, double *f_x, double *f_y,
                                  double *p, int PRINT_POLAR, double *IQUV) {
    int path_counter;
    double dl_current;

    double X_u[4], k_u[4], k_d[4];

    double Iinv, Iinv_pol;
    int POLARIZATION_ACTIVE = 0;

    double tetrad_u[4][4], tetrad_d[4][4];
    LOOP_ij tetrad_u[i][j] = 0.;
    LOOP_ij tetrad_d[i][j] = 0.;

    double photon_u_current[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
    double complex f_tetrad_u[4] = {0., 0., 0., 0.};
    double complex f_u[4] = {0., 0., 0., 0.};
    double complex S_A[4] = {0., 0., 0., 0.};

    struct GRMHD modvar;
    modvar.B = 0;
    modvar.n_e = 0.;
    modvar.theta_e = 0;

    LOOP_i {
        modvar.B_u[i] = 0;
        modvar.U_u[i] = 0;
        modvar.B_d[i] = 0;
        modvar.U_d[i] = 0;
    }
    modvar.igrid_c = -1;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        // PLASMA INTEGRATION STEP
        //////////////////////////

        double r_current = get_r(X_u);

        // Check whether the ray is currently in the GRMHD simulation volume
        if (get_fluid_params(X_u, &modvar) && r_current < OUTER_BOUND_POL) {
            pol_integration_step(modvar, frequency, &dl_current, C_CONST, X_u,
                                 k_u, k_d, &POLARIZATION_ACTIVE, f_u,
                                 f_tetrad_u, tetrad_d, tetrad_u, S_A, &Iinv,
                                 &Iinv_pol);
        } // End of if(IN_VOLUME)

        // SPACETIME-INTEGRATION STEP
        /////////////////////////////

        // If we HAVE been in-volume before, transport f_u (which is now
        // defined) one step. The final time this is done will be when
        // path_counter = 1; dl_current will then be at index 0 (path_counter -
        // 1).
        if (POLARIZATION_ACTIVE && path_counter > 0) {
            // Obtain the right k-vector, pointing back to observer, and
            // associated position. Pop into photon_u_current.
            LOOP_i {
                photon_u_current[i] = X_u[i];
                photon_u_current[i + 4] = k_u[i];
            }

            // One step: parallel transport of polarization vector.
            rk4_step_f(photon_u_current, f_u, dl_current);
        }
    } // End of for(path_counter...

    // CONSTRUCT FINAL (NON-INVARIANT) STOKES PARAMS SEEN BY OBSERVER
    /////////////////////////////////////////////////////////////////

    // Construct the observer tetrad.
    // X_u_current and k_u_current are simply the initial position and wave
    // vector. Note that k_u_current points INTO the camera sensor plane.
    LOOP_i {
        X_u[i] = lightpath[i];
        k_u[i] = lightpath[4 + i];
    }

    double complex f_obs_tetrad_u[4] = {0., 0., 0., 0.};
    construct_f_obs_tetrad_u(X_u, k_u, f_u, f_obs_tetrad_u);

    LOOP_i IQUV[i] = 0.;

    if (POLARIZATION_ACTIVE) {
        f_tetrad_to_stokes(Iinv, Iinv_pol, f_obs_tetrad_u, S_A);

        // Construct final (NON-INVARIANT) Stokes params.
        LOOP_i IQUV[i] = S_A[i] * pow(frequency, 3.);
    }
}
