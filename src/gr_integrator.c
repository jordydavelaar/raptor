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

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk45_step(double *y, void (*f)(double *, double *), double *dt, int bl) {
    // Array containing all "update elements" (4 times Nelements because RK4)
    double dx[DIM * 2 * 6];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK4
    double yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[DIM * 2];
    double err[DIM * 2];

    // Compute the RK4 update coefficients ('K_n' in lit., 'dx' here)
    int i, q;
    double weights[6] = {1. / 4., 3. / 8., 12. / 13.,
                         1.,      1. / 2., 0.}; // Weights used for updating y

    // slightly more cumbersome than rk4, since we now have non-zero crossterms
    // in the Butcher tablue
    // also U_u and A_u are shifted differently.
    q = 0;
    double b21 = 1. / 5.;

    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i];     // Use fvector to fill dx
        yshift[i] = y[i] + b21 * dx[0 * DIM * 2 + i]; // Update y
    }

    q = 1;
    double b31 = 3. / 40.;
    double b32 = 9. / 40.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b31 * dx[0 * DIM * 2 + i] +
                    b32 * dx[1 * DIM * 2 + i]; // Update y
    }

    q = 2;
    double b41 = 3. / 10.;
    double b42 = -9. / 10.;
    double b43 = 6. / 5.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b41 * dx[0 * DIM * 2 + i] +
                    b42 * dx[1 * DIM * 2 + i] +
                    b43 * dx[2 * DIM * 2 + i]; // Update y
    }

    q = 3;
    double b51 = -11. / 54.;
    double b52 = 5. / 2.;
    double b53 = -70. / 27.;
    double b54 = 35. / 27.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b51 * dx[0 * DIM * 2 + i] +
                    b52 * dx[1 * DIM * 2 + i] + b53 * dx[2 * DIM * 2 + i] +
                    b54 * dx[3 * DIM * 2 + i]; // Update y
    }

    q = 4;
    double b61 = 1631. / 55296.;
    double b62 = 175. / 512.;
    double b63 = 575. / 13824.;
    double b64 = 44275. / 110592.;
    double b65 = 253. / 4096.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b61 * dx[0 * DIM * 2 + i] +
                    b62 * dx[1 * DIM * 2 + i] + b63 * dx[2 * DIM * 2 + i] +
                    b64 * dx[3 * DIM * 2 + i] +
                    b65 * dx[4 * DIM * 2 + i]; // Update y
    }

    q = 5;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
    }

    double ch[6] = {37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.};
    double ct[6] = {2825. / 27648., 0.,     18575. / 48384., 13525. / 55296.,
                    277. / 14336.,  1. / 4.};

    // Update the y-vector (light ray)
    double errmax = -1;
    for (i = 0; i < DIM * 2; i++) {
        err[i] = fabs((ch[0] - ct[0]) * dx[0 * DIM * 2 + i] +
                      (ch[1] - ct[1]) * dx[1 * DIM * 2 + i] +
                      (ch[2] - ct[2]) * dx[2 * DIM * 2 + i] +
                      (ch[3] - ct[3]) * dx[3 * DIM * 2 + i] +
                      (ch[4] - ct[4]) * dx[4 * DIM * 2 + i] +
                      (ch[5] - ct[5]) * dx[5 * DIM * 2 + i]);
        if (err[i] > errmax)
            errmax = err[i];
    }

    double tol = 1e-3;

    if (errmax > tol && bl) {
        *dt = 0.9 * (*dt) * pow(tol / errmax, 1. / 5.);
        rk45_step(y, f, dt, 0);
    }

    if (errmax < tol / 10. && bl) {
        *dt = 1.5 * (*dt);
        rk45_step(y, f, dt, 0);
    }

    for (i = 0; i < DIM * 2; i++) {
        y[i] =
            y[i] + (ch[0] * dx[0 * DIM * 2 + i] + ch[1] * dx[1 * DIM * 2 + i] +
                    ch[2] * dx[2 * DIM * 2 + i] + ch[3] * dx[3 * DIM * 2 + i] +
                    ch[4] * dx[4 * DIM * 2 + i] + ch[5] * dx[5 * DIM * 2 + i]);
    }
}

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk4_step(double *y, void (*f)(double *, double *), double dt) {
    // Array containing all "update elements" (4 times Nelements because
    // RK4)
    double dx[DIM * 2 * 4];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK4
    double yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[DIM * 2];

    // Compute the RK4 update coefficients ('K_n' in lit., 'dx' here)
    int i, q;
    double weights[4] = {0.5, 0.5, 1., 0.}; // Weights used for updating y
    for (q = 0; q < 4; q++) {
        f(yshift,
          fvector); // Apply function f to current y to obtain fvector
        for (i = 0; i < DIM * 2; i++) {
            dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to fill dx
            yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < DIM * 2; i++)
        y[i] = y[i] + 1. / 6. *
                          (dx[0 * DIM * 2 + i] + dx[1 * DIM * 2 + i] * 2. +
                           dx[2 * DIM * 2 + i] * 2. + dx[3 * DIM * 2 + i]);
}

// Updates the vector y (containing position/velocity) by one RK2 step.
void rk2_step(double *y, void (*f)(double *, double *), double dt) {
    // Array containing all "update elements" (2 times Nelements because
    // RK2)
    double dx[DIM * 2 * 2];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK2
    double yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[DIM * 2];

    // Compute the RK2 update coefficients ('K_n' in lit., 'dx' here)
    int i, q;
    double weights[2] = {0.5, 0.}; // Weights used for updating y
    for (q = 0; q < 2; q++) {
        f(yshift,
          fvector); // Apply function f to current y to obtain fvector
        for (i = 0; i < DIM * 2; i++) {
            dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to update dx
            yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < DIM * 2; i++)
        y[i] = y[i] + dx[1 * DIM * 2 + i]; // y_n+1 = y_n + k2 + O(h^3)
}

// Updates the vector y (containing position/velocity) using 'velocity
// Verlet' Ref: Dolence et al 2009 eqn 14a-14d
void verlet_step(double *y, void (*f)(double *, double *), double dl) {
    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK2
    double yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[DIM * 2];

    // Temporary acceleration vector
    double A_u_temp[DIM];

    // Step 1: Compute A_u(lambda) (Preparation for Eq 14a)
    f(yshift, fvector); // fvector now contains A_u(lambda)

    // Step 2: Compute X_u(lambda + dlambda) and the temporary four-velocity
    // (Eq 14a, 14b)
    LOOP_i {
        yshift[i] += dl * yshift[i + DIM] + 0.5 * dl * dl * fvector[i + DIM];
        yshift[i + DIM] = yshift[i + DIM] + fvector[i + DIM] * dl;
        A_u_temp[i] = fvector[i + DIM]; // STORE A_u(lambda)
    }

    // Step 3: Compute A_u(lambda + dlambda) (Eq 14c)
    f(yshift, fvector); // fvector now contains A_u(lambda + dl)

    // Step 4: Compute new velocity (Eq 14d)
    LOOP_i {
        y[i] = yshift[i]; // X_u(l + dl)
        y[i + DIM] += 0.5 * (A_u_temp[i] + fvector[i + DIM]) * dl; // A_u(l+dl)
    }
}

// Returns an appropriate stepsize dlambda, which depends on position &
// velocity Ref. DOLENCE & MOSCIBRODZKA 2009
double stepsize(double X_u[4], double U_u[4]) {
    double SMALL = 1.e-40;
#if (metric == CKS)
    double dlx1 = STEPSIZE / (fabs(U_u[1]) + SMALL * SMALL);
    double dlx2 = STEPSIZE / (fabs(U_u[2]) + SMALL * SMALL);
    double dlx3 = STEPSIZE / (fabs(U_u[3]) + SMALL * SMALL);
#else
#if (metric == MKSHARM)
    double pf = 1.;
#else
    double pf = M_PI;
#endif

    double dlx1 = STEPSIZE / (fabs(U_u[1]) + SMALL * SMALL);
    double dlx2 =
        STEPSIZE * fmin(X_u[2], pf - X_u[2]) / (fabs(U_u[2]) + SMALL * SMALL);
    double dlx3 = STEPSIZE / (fabs(U_u[3]) + SMALL * SMALL);
#endif
    double idlx1 = 1. / (fabs(dlx1) + SMALL * SMALL);
    double idlx2 = 1. / (fabs(dlx2) + SMALL * SMALL);
    double idlx3 = 1. / (fabs(dlx3) + SMALL * SMALL);
#if (metric == CKS)
    double r = get_r(X_u);
    return -sqrt(r) / (idlx1 + idlx2 + idlx3);
#else
    return -1. / (idlx1 + idlx2 + idlx3);
#endif
}

// The function to be used by the integrator for GR geodesic calculations.
// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_geodesic(double *y, double *fvector) {
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];

    // Initialize position, four-velocity, and four-acceleration vectors
    // based on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double A_u[4] = {0., 0., 0., 0.};         // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
    // connection_udd(X_u, gamma_udd);
    connection_num_udd(X_u, gamma_udd);

    // Compute 4-acceleration using the geodesic equation
    // Einstein summation over indices v and w
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];

    // Update fvector
    LOOP_i {
        fvector[i] = U_u[i];
        fvector[i + DIM] = A_u[i];
    }
}

// Integrate the null geodesic defined by "photon_u"
void integrate_geodesic(double alpha, double beta, double *lightpath,
                        int *steps, double cutoff_inner) {
    int q;
    double t_init = 0.;
    double dlambda_adaptive = -1.0;
    int theta_turns = 0;
    double thetadot_prev;
    double X_u[4], k_u[4];
    double photon_u[8];
    // Create initial ray conditions
    initialize_photon(alpha, beta, photon_u, t_init);
    LOOP_i X_u[i] = photon_u[i];
    // Current r-coordinate
    double r_current = get_r(X_u);

    // Reset lambda and steps
    double lambda = 0.;
    *steps = 0;

    int TERMINATE = 0; // Termination condition for ray

    // Trace light ray until it reaches the event horizon or the outer
    // cutoff, or steps > max_steps
#if (metric == BL || metric == MBL)

    // Stop condition for BL coords
    while (r_current > cutoff_inner && r_current < cutoff_outer &&
           *steps < max_steps && !TERMINATE) { // && photon_u[0] < t_final){

#else

    // Stop condition for KS coords
    while (r_current < cutoff_outer && r_current > cutoff_inner &&
           *steps < max_steps && !TERMINATE) {

#endif

        // Current photon position/wave vector
        LOOP_i {
            X_u[i] = photon_u[i];
            k_u[i] = photon_u[i + 4];
        }

        // Possibly terminate ray to eliminate higher order images
        if (thetadot_prev * photon_u[6] < 0. && *steps > 2)
            theta_turns += 1;
        thetadot_prev = photon_u[6];
        if ((beta < 0. && theta_turns > max_order) ||
            (beta > 0. && theta_turns > (max_order + 1)))
            TERMINATE = 1;

// Compute adaptive step size
// dlambda_adaptive = -STEPSIZE;
#if (int_method != RK45)
        dlambda_adaptive = stepsize(X_u, k_u);
#endif
        // Enter current position/velocity/dlambda into lightpath
        for (q = 0; q < 8; q++)
            lightpath[*steps * 9 + q] = photon_u[q];
        lightpath[*steps * 9 + 8] = fabs(dlambda_adaptive);

        // Advance ray/particle
#if (int_method == RK4)

        rk4_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == VER)

        verlet_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == RK45)

    rk45_step(photon_u, &f_geodesic, &dlambda_adaptive, 1);
#endif
        LOOP_i X_u[i] = photon_u[i];
        LOOP_i k_u[i] = photon_u[i + 4];
        // Advance (affine) parameter lambda
        lambda += fabs(dlambda_adaptive);
        r_current = get_r(X_u);

        *steps = *steps + 1;
    }
}
