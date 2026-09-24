/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <gsl/gsl_sf_ellint.h>

// FUNCTIONS
////////////

// Updates y with a recursive RKF method. Combines a 4th and 5th order RK step
// Uses error estimate to perform an updated step if error too large.
void rk45_step(double *y, void (*f)(double *, double *), double *dt, int bl) {
    // Array containing all "update elements" (4 times Nelements because RK4)
    double dx[DIM * 2 * 6];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK4
    double yshift[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};
    double yold[DIM * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[DIM * 2];
    double err[DIM * 2];

    // slightly more cumbersome than rk4, since we now have non-zero crossterms
    // in the Butcher tablue
    // also U_u and A_u are shifted differently.
    int q = 0;
    double b21 = 1. / 5.;

    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (int i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i];     // Use fvector to fill dx
        yshift[i] = y[i] + b21 * dx[0 * DIM * 2 + i]; // Update y
    }

    q = 1;
    double b31 = 3. / 40.;
    double b32 = 9. / 40.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (int i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b31 * dx[0 * DIM * 2 + i] +
                    b32 * dx[1 * DIM * 2 + i]; // Update y
    }

    q = 2;
    double b41 = 3. / 10.;
    double b42 = -9. / 10.;
    double b43 = 6. / 5.;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (int i = 0; i < DIM * 2; i++) {
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
    for (int i = 0; i < DIM * 2; i++) {
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
    for (int i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
        yshift[i] = y[i] + b61 * dx[0 * DIM * 2 + i] +
                    b62 * dx[1 * DIM * 2 + i] + b63 * dx[2 * DIM * 2 + i] +
                    b64 * dx[3 * DIM * 2 + i] +
                    b65 * dx[4 * DIM * 2 + i]; // Update y
    }

    q = 5;
    f(yshift, fvector); // Apply function f to current y to obtain fvector
    for (int i = 0; i < DIM * 2; i++) {
        dx[q * DIM * 2 + i] = (*dt) * fvector[i]; // Use fvector to fill dx
    }

    double ch[6] = {37. / 378., 0, 250. / 621., 125. / 594., 0., 512. / 1771.};
    double ct[6] = {2825. / 27648., 0.,     18575. / 48384., 13525. / 55296.,
                    277. / 14336.,  1. / 4.};

    // Update the y-vector (light ray)

    double tol = 1e-9;

    double errmax = -1;
    for (int i = 0; i < DIM * 2; i++) {
        err[i] = fabs((ch[0] - ct[0]) * dx[0 * DIM * 2 + i] +
                      (ch[1] - ct[1]) * dx[1 * DIM * 2 + i] +
                      (ch[2] - ct[2]) * dx[2 * DIM * 2 + i] +
                      (ch[3] - ct[3]) * dx[3 * DIM * 2 + i] +
                      (ch[4] - ct[4]) * dx[4 * DIM * 2 + i] +
                      (ch[5] - ct[5]) * dx[5 * DIM * 2 + i]);
    }

    for (int i = 0; i < DIM * 2; i++) {
        y[i] =
            y[i] + (ch[0] * dx[0 * DIM * 2 + i] + ch[1] * dx[1 * DIM * 2 + i] +
                    ch[2] * dx[2 * DIM * 2 + i] + ch[3] * dx[3 * DIM * 2 + i] +
                    ch[4] * dx[4 * DIM * 2 + i] + ch[5] * dx[5 * DIM * 2 + i]);
    }

    for (int i = 0; i < DIM; i++) {
        double yscal =
            fabs(ch[0] * dx[0 * DIM * 2 + i] + ch[1] * dx[1 * DIM * 2 + i] +
                 ch[2] * dx[2 * DIM * 2 + i] + ch[3] * dx[3 * DIM * 2 + i] +
                 ch[4] * dx[4 * DIM * 2 + i] + ch[5] * dx[5 * DIM * 2 + i]);

        err[i] = err[i] / yscal;
        if (err[i] > errmax)
            errmax = err[i];
    }
    //    errmax/=tol;

    if (errmax > tol && bl) {
        (*dt) = -fmax(0.001 * fabs(*dt),
                      0.84 * fabs(*dt) * powf(tol / errmax, 1. / 4.));
        rk45_step(yold, f, dt, 0);

        for (int i = 0; i < DIM * 2; i++) {
            y[i] = yold[i];
        }
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
    double SMALL = 1.e-80;
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
    return -(r) / (idlx1 + idlx2 + idlx3);
#else

#if(metric==MKSBHAC)

       double step_grid =1e100;
       if(exp(X_u[1])<RT_OUTER_CUTOFF){
           int igrid = find_igrid(X_u, block_info, Xgrid);
	   double step1 = block_info[igrid].dxc_block[0]/(fabs(U_u[1]) + SMALL * SMALL);
           double step2 = block_info[igrid].dxc_block[1]/(fabs(U_u[2]) + SMALL * SMALL);
           double step3 = block_info[igrid].dxc_block[2]/(fabs(U_u[3]) + SMALL * SMALL);
           double istep1 = 1. / (fabs(step1) + SMALL * SMALL);
           double istep2 = 1. / (fabs(step2) + SMALL * SMALL);
           double istep3 = 1. / (fabs(step3) + SMALL * SMALL);
           double minstep = 1/(istep1+istep2+istep3);

           step_grid =  minstep/4.;
        }
        double step = fmin(1/ (idlx1 + idlx2 + idlx3), step_grid);
        return -step;

#else
    return -1. / (idlx1 + idlx2 + idlx3);
#endif

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
#if (metric == MKSBHAC || metric == MKSHARM)
    connection_udd(X_u, gamma_udd);
#else
    connection_num_udd(X_u, gamma_udd);
#endif
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
    double dlambda_adaptive = -0.1;
    int theta_turns = 0;
    double thetadot_prev;
    double X_u[4], k_u[4];
    double photon_u[8];
    // Create initial ray conditions
    initialize_photon(alpha, beta, photon_u, t_init);
    LOOP_i X_u[i] = photon_u[i];
    LOOP_i k_u[i] = photon_u[i+4];
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

        LOOP_i {
            photon_u[i+4] = k_u[i];
        }

        // Enter current position/velocity/dlambda into lightpath
        for (q = 0; q < 8; q++)
            lightpath[*steps * 9 + q] = photon_u[q];

        // Possibly terminate ray to eliminate higher order images
        if (thetadot_prev * photon_u[6] < 0. && *steps > 2)
            theta_turns += 1;
        thetadot_prev = photon_u[6];
        if ((beta < 0. && theta_turns > max_order) ||
            (beta > 0. && theta_turns > (max_order + 1)))
            TERMINATE = 1;

        // Compute educational guess for an adaptive step size
        // dlambda_adaptive = -STEPSIZE;
        dlambda_adaptive = stepsize(X_u, k_u);

        // Advance ray/particle
#if (int_method == RK2)

        rk2_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == RK4)

        rk4_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == VER)

    verlet_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == RK45)

    rk45_step(photon_u, &f_geodesic, &dlambda_adaptive, 1);
#endif


        lightpath[*steps * 9 + 8] = fabs(dlambda_adaptive);

        LOOP_i X_u[i] = photon_u[i];
        LOOP_i k_u[i] = photon_u[i + 4];

        // Advance (affine) parameter lambda
        lambda += fabs(dlambda_adaptive);
        r_current = get_r(X_u);

        *steps = *steps + 1;
    }
}

// Winding number n = Delta_phi / (2*pi) accumulated along the geodesic.
// phi (lightpath index 3 within each 9-entry step) is stored unwrapped by
// integrate_geodesic, so the total signed sweep is just the difference
// between the last and first stored phi values.
double compute_photon_order(double *lightpath, int steps) {
    if (steps < 1)
        return 0.;
    double dphi_total = lightpath[(steps - 1) * 9 + 3] - lightpath[0 * 9 + 3];
    return dphi_total / (2. * M_PI);
}

// Complete elliptic integral of the first kind, K(m), parameter convention
// m = k^2 (matches Python/scipy's ellipk(m)), valid for m<1 including
// negative m via the imaginary-modulus transform K(m) = K(m/(m-1))/sqrt(1-m)
// (verified numerically against scipy.special.ellipk). GSL's
// gsl_sf_ellint_Kcomp takes the modulus k, not m, hence the sqrt()s below.
static double ellipK_param(double m) {
    if (m >= 1.)
        m = 1. - 1e-12;
    if (m < 0.) {
        double mt = m / (m - 1.);
        return gsl_sf_ellint_Kcomp(sqrt(mt), GSL_PREC_DOUBLE) / sqrt(1. - m);
    }
    return gsl_sf_ellint_Kcomp(sqrt(m), GSL_PREC_DOUBLE);
}

// Two additional photon-subring diagnostics computed from the already
// integrated geodesic:
//
// *n_eq_crossings: number of times the ray crosses the equatorial plane
// (theta = pi/2). Simple integer order label used e.g. by Gralla, Lupsasca &
// Marrone 2020 (PRD 102, 124004) and by adaptive-ray-tracing lensing-band
// definitions (m=0 direct image, m=1 first indirect, ...).
//
// *mino_order: continuous half-orbit count following the photon's polar
// (theta) motion, n = tau_mino * sqrt(-a^2 u_minus) / (2*K(u_plus/u_minus)),
// where tau_mino = integral d(affine parameter)/Sigma is the elapsed Mino
// time and u_plus/u_minus follow Gralla & Lupsasca 2020 (arXiv:1910.12873,
// hereafter GL20), Eq. 11.
//
// IMPORTANT: the "2*K" here is NOT the "4*K" that appears in GL20 Eq. 36
// itself (also implemented as kgeo's n_poloidal_orbits). GL20 Eq. 36 defines
// n_GL36 = tau_mino*sqrt(-a^2 u_minus)/(4*K), the number of FULL poloidal
// orbits (one full theta_min -> theta_max -> theta_min cycle increments it
// by 1). What this function returns is 2*n_GL36, i.e. the number of HALF
// orbits, because that is the quantity conventionally used to label photon
// subrings n=1,2,3,... (Delta(mino_order)=1 per turning-point-to-turning-
// point interval). This bridges two separately-stated facts rather than
// being read off a single equation: GL20 Eq. 36 for n_GL36 itself, plus
// Himwich, Johnson, Lupsasca & Strominger 2020 (arXiv:2001.08750), which
// states in its introduction that "the n-th subring is comprised of photons
// that circumnavigate the black hole n/2 times" -- i.e. n_subring =
// 2 * (number of full orbits) = 2 * n_GL36, which is exactly the factor of 2
// (not 4) used below. To recover the literal GL20 Eq. 36 quantity instead,
// divide this function's mino_order output by 2.
//
// Undefined (sentinel -999) for vortical geodesics (eta<=0), which never
// cross the equatorial plane and have no polar libration to count fractions
// of -- this includes every ray for a camera placed exactly edge-on
// (inclination 90 deg) at zero impact parameter beta.
//
// lambda, eta reuse the exact (E=1) conserved-quantity expressions from
// initialize_photon() in metric.c, for self-consistency with whatever
// geodesic RAPTOR actually integrates for this pixel.
void compute_photon_order_extra(double *lightpath, int steps, double alpha,
                                double beta, double *mino_order,
                                int *n_eq_crossings) {
    *n_eq_crossings = 0;
    *mino_order = -999.;

    if (steps < 2)
        return;

#if (metric == MKSBHAC || metric == MKSHARM)
    double eq_ref = 0.5; // X2 = 0.5 <=> theta = pi/2 for this coordinate map
#else
    double eq_ref = M_PI / 2.;
#endif

    double tau_mino = 0.;
    double prev = lightpath[0 * 9 + 2] - eq_ref;
    int crossings = 0;
    for (int q = 0; q < steps - 1; q++) {
        double r_q = logscale ? exp(lightpath[q * 9 + 1]) : lightpath[q * 9 + 1];
        double th2_q = lightpath[q * 9 + 2];
#if (metric == MKSBHAC || metric == MKSHARM)
        double theta_q =
            M_PI * th2_q + 0.5 * (1. - hslope) * sin(2. * M_PI * th2_q);
#else
        double theta_q = th2_q;
#endif
        double sigma_q = r_q * r_q + a * a * cos(theta_q) * cos(theta_q);
        double dlambda_q = lightpath[q * 9 + 8];
        tau_mino += dlambda_q / sigma_q;

        double cur = lightpath[(q + 1) * 9 + 2] - eq_ref;
        if (prev * cur < 0.)
            crossings++;
        prev = cur;
    }
    *n_eq_crossings = crossings;

    double mu0 = cos(INCLINATION / 180. * M_PI);
    double lam = -alpha * sqrt(1. - mu0 * mu0);
    double eta = beta * beta + mu0 * mu0 * (alpha * alpha - 1.);

    if (eta <= 0.)
        return; // vortical geodesic: no polar libration to normalize by

    double u_plus, u_minus, a2u_minus;
    if (fabs(a) < 1e-6) {
        u_plus = eta / (eta + lam * lam);
        u_minus = u_plus;
        a2u_minus = -(eta + lam * lam);
    } else {
        double a2 = a * a;
        double Delta_theta = 0.5 * (1. - (eta + lam * lam) / a2);
        double disc = Delta_theta * Delta_theta + eta / a2;
        disc = disc > 0. ? sqrt(disc) : 0.;
        u_plus = Delta_theta + disc;
        u_minus = Delta_theta - disc;
        a2u_minus = a2 * u_minus;
    }

    if (a2u_minus >= 0. || u_minus == 0.)
        return; // degenerate: no polar libration to normalize by

    double uratio = u_plus / u_minus;
    double K = ellipK_param(uratio);
    double half_period = 2. * K / sqrt(-a2u_minus);

    *mino_order = tau_mino / half_period;
}
