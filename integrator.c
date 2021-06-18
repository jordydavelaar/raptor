/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 */

#include "constants.h"
#include "functions.h"
#include "parameters.h"
#include <complex.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk4_step(double *y, void (*f)(double *, double *), double dt) {
    // Array containing all "update elements" (4 times Nelements because RK4)
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
        f(yshift, fvector); // Apply function f to current y to obtain fvector
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
    // Array containing all "update elements" (2 times Nelements because RK2)
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
        f(yshift, fvector); // Apply function f to current y to obtain fvector
        for (i = 0; i < DIM * 2; i++) {
            dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to update dx
            yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < DIM * 2; i++)
        y[i] = y[i] + dx[1 * DIM * 2 + i]; // y_n+1 = y_n + k2 + O(h^3)
}

// Updates the vector y (containing position/velocity) using 'velocity Verlet'
// Ref: Dolence et al 2009 eqn 14a-14d
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
    int i;
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

// Returns an appropriate stepsize dlambda, which depends on position & velocity
// Ref. DOLENCE & MOSCIBRODZKA 2009
double stepsize(double X_u[4], double U_u[4]) {
    double SMALL = 1.e-40;

    double dlx1 = STEPSIZE / (fabs(U_u[1]) + SMALL * SMALL);
    double dlx2 =
        STEPSIZE * fmin(X_u[2], M_PI - X_u[2]) / (fabs(U_u[2]) + SMALL * SMALL);
    double dlx3 = STEPSIZE / (fabs(U_u[3]) + SMALL * SMALL);

    double idlx1 = 1. / (fabs(dlx1) + SMALL * SMALL);
    double idlx2 = 1. / (fabs(dlx2) + SMALL * SMALL);
    double idlx3 = 1. / (fabs(dlx3) + SMALL * SMALL);

    return -1. / (idlx1 + idlx2 + idlx3);
}

// The function to be used by the integrator for GR geodesic calculations.
// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_geodesic(double *y, double *fvector) {
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double A_u[4] = {0., 0., 0., 0.};         // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
    connection_udd(X_u, gamma_udd);
    // connection_num_udd(X_u, gamma_udd);

    // Compute 4-acceleration using the geodesic equation
    int i, j, k; // Einstein summation over indices v and w
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];

    // Update fvector
    LOOP_i {
        fvector[i] = U_u[i];
        fvector[i + DIM] = A_u[i];
    }
}

// Integrate the null geodesic defined by "photon_u"
void integrate_geodesic(double alpha, double beta, double *photon_u,
                        double *lightpath, int *steps, double cutoff_inner) {
    int i, q;
    double t_init = 0.;
    double dlambda_adaptive;
    int theta_turns = 0;
    double thetadot_prev;
    double X_u[4], k_u[4];

    // Create initial ray conditions
    initialize_photon(alpha, beta, photon_u, t_init);

    // Current r-coordinate
    double r_current = logscale ? exp(photon_u[1]) : photon_u[1];

    // Reset lambda and steps
    double lambda = 0.;
    *steps = 0;

    int TERMINATE = 0; // Termination condition for ray

    // Trace light ray until it reaches the event horizon or the outer
    // cutoff, or steps > max_steps
#if (metric == BL || metric == MBL || metric == DM)

    // Stop condition for BL coords
    while (r_current > cutoff_inner && r_current < cutoff_outer &&
           *steps < max_steps && !TERMINATE) { // && photon_u[0] < t_final){

#elif (metric == KS || metric == MKS || metric == MKS2)

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
        dlambda_adaptive = stepsize(X_u, k_u);

        // Enter current position/velocity/dlambda into lightpath
        for (q = 0; q < 8; q++)
            lightpath[*steps * 9 + q] = photon_u[q];
        lightpath[*steps * 9 + 8] = fabs(dlambda_adaptive);

        // Advance ray/particle
#if (int_method == RK4)

        rk4_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == VER)

        verlet_step(photon_u, &f_geodesic, dlambda_adaptive);

#endif

        // Advance (affine) parameter lambda
        lambda += fabs(dlambda_adaptive);
        r_current = logscale ? exp(photon_u[1]) : photon_u[1];

        *steps = *steps + 1;
    }
}


/*
// Recursively calculates the determinant of a matrix.
// Source: http://ideone.com/fork/92JF0O
double determ(double matrix[][4], int n) {
    int p, h, k, i, j;
    double det = 0.;
    double temp[4][4];
    if (n == 1) {
        return matrix[0][0];
    } else if (n == 2) {
        det = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        return det;
    } else {
        for (p = 0; p < n; p++) {
            h = 0;
            k = 0;
            for (i = 1; i < n; i++)
                for (j = 0; j < n; j++) {
                    if (j == p)
                        continue;
                    temp[h][k] = matrix[i][j];
                    k++;
                    if (k == n - 1) {
                        h++;
                        k = 0;
                    }
                }
            det = det + matrix[0][p] * pow(-1, p) * determ(temp, n - 1);
        }
        return det;
    }
}

// Here we implement Ziri's tetrad
void create_tetrad_u2(const double X_u[], const double k_u[],
                      const double U_u[], double tetrad_u[][4]) {
    // Summation indices:
    int i, j, k, l;

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    // IN GAMMIE 2012 NOTATION:
    double e_u_t[4];
    LOOP_i e_u_t[i] = U_u[i]; // Check

    double e_u_K[4];
    double omega = -inner_product(X_u, k_u, U_u); // Ziri calls omega "g"
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i];    // Check

    double e_u_para[4];

    double e_u_perp[4];
    double U_d[4];
    double k_d[4];

    LOOP_i U_d[i] = 0.;
    LOOP_i k_d[i] = 0.;

    // Need U_d, k_d:
    lower_index(X_u, U_u, U_d);
    lower_index(X_u, k_u, k_d);

    // TRIAL VECTOR b_u
    ///////////////////

    // Only requirement: b dot b > 0 (i.e., b is spacelike)

    double b_u[4] = {0., 0., 0., 0.};

    // Strategy for creating b:

    // 1. Pick random numbers between 0 and 1.
    LOOP_i b_u[i] = 1. / sqrt(3.) * genrand_RCARRY();

    // 2. Check inner product.
    // fprintf(stderr, "\nb dot b = %+.15e", inner_product(X_u, b_u, b_u));

    // 3. If spacelike, OK.
    // 4. If not, flip sign of b_u[0]? Or just create new b_u?
    while (inner_product(X_u, b_u, b_u) < 0.1) {
        LOOP_i b_u[i] = 1. / sqrt(3.) * genrand_RCARRY();
    }

    // Now we have b_u s.t. b dot b is spacelike. CHECK!

    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);

    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] =
        (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive; // CHECK

    // Permutation symbol eta, [ijkl]. Can be made into the contravariant
    // Levi-Civita tensor by multiplying with -1/sqrt(g), g being the
    // determinant of the covariant metric. (source:
    // cpluscode.blogspot.nl/p/physics-c.html)
    double g = determ(g_dd, 4);
    double eps[4][4][4][4];
    LOOP_ijkl {
        if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) ||
            (k == l))
            eps[i][j][k][l] = 0.;
        else
            eps[i][j][k][l] = ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) *
                               (k - l) / 12.);
    }

    // fprintf(stderr, "\n\nU DOT T: %+.15e\n\n", inner_product(X_u, U_u, t_u));
    // fprintf(stderr, "k DOT T: %+.15e\n\n", inner_product(X_u, k_u, t_u));
    // fprintf(stderr, "k DOT k: %+.15e\n\n", inner_product(X_u, k_u, k_u));
    // fprintf(stderr, "U DOT U: %+.15e\n\n", inner_product(X_u, U_u, U_u));

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] +=
        (-1. / sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) /
        (omega * Ncursive); // Check

    // Construct the tetrad with contravariant coordinate index
    // CONVENTION: t, para, perp, K <=> t, x, y, z
    tetrad_u[0][0] = e_u_t[0];
    tetrad_u[1][0] = e_u_t[1];
    tetrad_u[2][0] = e_u_t[2];
    tetrad_u[3][0] = e_u_t[3];

    tetrad_u[0][3] = e_u_K[0];
    tetrad_u[1][3] = e_u_K[1];
    tetrad_u[2][3] = e_u_K[2];
    tetrad_u[3][3] = e_u_K[3];

    tetrad_u[0][1] = e_u_para[0];
    tetrad_u[1][1] = e_u_para[1];
    tetrad_u[2][1] = e_u_para[2];
    tetrad_u[3][1] = e_u_para[3];

    tetrad_u[0][2] = e_u_perp[0];
    tetrad_u[1][2] = e_u_perp[1];
    tetrad_u[2][2] = e_u_perp[2];
    tetrad_u[3][2] = e_u_perp[3];
}

// Here we implement Ziri's tetrad
void create_observer_tetrad_u2(const double X_u[], const double k_u[],
                               const double U_u[], const double b_u[],
                               double tetrad_u[][4]) {
    // Summation indices:
    int i, j, k, l;

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    // IN GAMMIE 2012 NOTATION:
    double e_u_t[4];
    LOOP_i e_u_t[i] = U_u[i]; // Check

    double e_u_K[4];
    double omega = -inner_product(X_u, k_u, U_u); // Ziri calls omega "g"
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i];    // Check

    double e_u_para[4];

    double e_u_perp[4];
    double U_d[4];
    double k_d[4];

    LOOP_i U_d[i] = 0.;
    LOOP_i k_d[i] = 0.;

    // Need U_d, k_d:
    lower_index(X_u, U_u, U_d);
    lower_index(X_u, k_u, k_d);

    // CAM-UP VECTOR
    ////////////////

    // Only requirement: b dot b > 0 (i.e., b is spacelike) and b points along Y
    // direction of image.

    //    double b_u[4] = {0., 0., 1., 0.};

    // Now we have b_u s.t. b dot b is spacelike. CHECK!

    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);

    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] =
        (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive; // CHECK

    // Permutation symbol eta, [ijkl]. Can be made into the contravariant
    // Levi-Civita tensor by multiplying with -1/sqrt(g), g being the
    // determinant of the covariant metric. (source:
    // cpluscode.blogspot.nl/p/physics-c.html)
    double g = determ(g_dd, 4);
    double eps[4][4][4][4];
    LOOP_ijkl {
        if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) ||
            (k == l))
            eps[i][j][k][l] = 0.;
        else
            eps[i][j][k][l] = ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) *
                               (k - l) / 12.);
    }

    // fprintf(stderr, "\n\nU DOT T: %+.15e\n\n", inner_product(X_u, U_u, t_u));
    // fprintf(stderr, "k DOT T: %+.15e\n\n", inner_product(X_u, k_u, t_u));
    // fprintf(stderr, "k DOT k: %+.15e\n\n", inner_product(X_u, k_u, k_u));
    // fprintf(stderr, "U DOT U: %+.15e\n\n", inner_product(X_u, U_u, U_u));

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] +=
        (-1. / sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) /
        (omega * Ncursive); // Check

    // Construct the tetrad with contravariant coordinate index
    // CONVENTION: t, para, perp, K <=> t, x, y, z
    tetrad_u[0][0] = e_u_t[0];
    tetrad_u[1][0] = e_u_t[1];
    tetrad_u[2][0] = e_u_t[2];
    tetrad_u[3][0] = e_u_t[3];

    tetrad_u[0][3] = e_u_K[0];
    tetrad_u[1][3] = e_u_K[1];
    tetrad_u[2][3] = e_u_K[2];
    tetrad_u[3][3] = e_u_K[3];

    tetrad_u[0][2] = e_u_para[0];
    tetrad_u[1][2] = e_u_para[1];
    tetrad_u[2][2] = e_u_para[2];
    tetrad_u[3][2] = e_u_para[3];

    tetrad_u[0][1] = e_u_perp[0];
    tetrad_u[1][1] = e_u_perp[1];
    tetrad_u[2][1] = e_u_perp[2];
    tetrad_u[3][1] = e_u_perp[3];
}

double tetrad_identity_eta(const double X_u[4], const double tetrad_u[4][4],
                           const int a, const int b) {
    double result = 0.;

    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    int i, j;
    LOOP_ij result += g_dd[i][j] * tetrad_u[i][a] * tetrad_u[j][b];

    return result;
}

double tetrad_identity_g(const double tetrad_u[][4], const int mu,
                         const int nu) {

    double eta_uu[4][4] = {{-1., 0., 0., 0.},
                           {0., 1., 0., 0.},
                           {0., 0., 1., 0.},
                           {0., 0., 0., 1.}};

    double result = 0.;

    int i, j;
    LOOP_ij result += eta_uu[i][j] * tetrad_u[mu][i] * tetrad_u[nu][j];

    return result;
}

double tetrad_identity_sum_latin(const double tetrad_u[4][4],
                                 const double tetrad_d[4][4], const int mu,
                                 const int nu) {

    double result = 0.;

    int i;
    LOOP_i result += tetrad_u[mu][i] * tetrad_d[nu][i];

    return result;
}

double tetrad_identity_sum_greek(const double tetrad_u[4][4],
                                 const double tetrad_d[4][4], const int a,
                                 const int b) {

    double result = 0.;

    int i;
    LOOP_i result += tetrad_u[i][a] * tetrad_d[i][b];

    return result;
}

void create_tetrad_d(const double X_u[], const double tetrad_u[][4],
                     double tetrad_d[][4]) {
    double eta_minkowski[4][4] = {
        {-1., 0., 0., 0.},
        {0., 1., 0., 0.},
        {0., 0., 1., 0.},
        {0., 0., 0., 1.},
    };

    int i, j, q, s, k, l;

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    // Create the tetrad with covariant coordinate index:
    LOOP_ij tetrad_d[i][j] = 0.;

    // ***NOTE*** the index order must be swapped here, i.e. we take the
    // transpose.
    LOOP_ij {
        LOOP_kl tetrad_d[j][i] +=
            eta_minkowski[i][k] * g_dd[j][l] * tetrad_u[l][k];
    }
}

double check_tetrad_compact(const double X_u[], const double tetrad_u[][4]) {
    double result = 0.;

    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 0) - (-1.)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 1)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 2)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 0)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 1) - 1.) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 2)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 0)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 1)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 2) - 1.) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 0)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 1)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 2)) +
              fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 3) - 1.);

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double tetrad_d[4][4];
    int i, j;
    LOOP_ij tetrad_d[i][j] = 0.;
    create_tetrad_d(X_u, tetrad_u, tetrad_d);

    result += fabs(tetrad_identity_g(tetrad_u, 0, 0) - g_uu[0][0]) +
              fabs(tetrad_identity_g(tetrad_u, 0, 1) - g_uu[0][1]) +
              fabs(tetrad_identity_g(tetrad_u, 0, 2) - g_uu[0][2]) +
              fabs(tetrad_identity_g(tetrad_u, 0, 3) - g_uu[0][3]);
    result += fabs(tetrad_identity_g(tetrad_u, 1, 0) - g_uu[1][0]) +
              fabs(tetrad_identity_g(tetrad_u, 1, 1) - g_uu[1][1]) +
              fabs(tetrad_identity_g(tetrad_u, 1, 2) - g_uu[1][2]) +
              fabs(tetrad_identity_g(tetrad_u, 1, 3) - g_uu[1][3]);
    result += fabs(tetrad_identity_g(tetrad_u, 2, 0) - g_uu[2][0]) +
              fabs(tetrad_identity_g(tetrad_u, 2, 1) - g_uu[2][1]) +
              fabs(tetrad_identity_g(tetrad_u, 2, 2) - g_uu[2][2]) +
              fabs(tetrad_identity_g(tetrad_u, 2, 3) - g_uu[2][3]);
    result += fabs(tetrad_identity_g(tetrad_u, 3, 0) - g_uu[3][0]) +
              fabs(tetrad_identity_g(tetrad_u, 3, 1) - g_uu[3][1]) +
              fabs(tetrad_identity_g(tetrad_u, 3, 2) - g_uu[3][2]) +
              fabs(tetrad_identity_g(tetrad_u, 3, 3) - g_uu[3][3]);

    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 0) - 1.) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 1)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 2)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 0)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 1) - 1.) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 2)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 0)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 1)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 2) - 1.) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 0)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 1)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 2)) +
              fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 3) - 1.);

    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 0) - 1.) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 1)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 2)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 0)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 1) - 1.) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 2)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 0)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 1)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 2) - 1.) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 0)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 1)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 2)) +
              fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 3) - 1.);
    if (isnan(result))
        fprintf(stderr, "position of nan is %e %e %e\n", exp(X_u[1]), X_u[2],
                X_u[3]);
    // printf("\nTetrad identities (should be close to zero): %+.15e", result);
    return result;
}

void check_tetrad_identities(const double X_u[], double tetrad_u[][4]) {
    printf("\nRecover dirac delta");

    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 0) - (-1.)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 1)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 2)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 3)));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 0)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 1) - 1.),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 2)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 3)));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 0)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 1)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 2) - 1.),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 3)));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 0)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 1)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 2)),
           fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 3) - 1.));

    printf("\nRecover metric");

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double tetrad_d[4][4];
    create_tetrad_d(X_u, tetrad_u, tetrad_d);

    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_g(tetrad_u, 0, 0) - g_uu[0][0]),
           fabs(tetrad_identity_g(tetrad_u, 0, 1) - g_uu[0][1]),
           fabs(tetrad_identity_g(tetrad_u, 0, 2) - g_uu[0][2]),
           fabs(tetrad_identity_g(tetrad_u, 0, 3) - g_uu[0][3]));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_g(tetrad_u, 1, 0) - g_uu[1][0]),
           fabs(tetrad_identity_g(tetrad_u, 1, 1) - g_uu[1][1]),
           fabs(tetrad_identity_g(tetrad_u, 1, 2) - g_uu[1][2]),
           fabs(tetrad_identity_g(tetrad_u, 1, 3) - g_uu[1][3]));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_g(tetrad_u, 2, 0) - g_uu[2][0]),
           fabs(tetrad_identity_g(tetrad_u, 2, 1) - g_uu[2][1]),
           fabs(tetrad_identity_g(tetrad_u, 2, 2) - g_uu[2][2]),
           fabs(tetrad_identity_g(tetrad_u, 2, 3) - g_uu[2][3]));
    printf("\n%.3e %.3e %.3e %.3e",
           fabs(tetrad_identity_g(tetrad_u, 3, 0) - g_uu[3][0]),
           fabs(tetrad_identity_g(tetrad_u, 3, 1) - g_uu[3][1]),
           fabs(tetrad_identity_g(tetrad_u, 3, 2) - g_uu[3][2]),
           fabs(tetrad_identity_g(tetrad_u, 3, 3) - g_uu[3][3]));

    printf("\n");

    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 0),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 1),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 2),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 0),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 1),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 2),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 0),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 1),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 2),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 0),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 1),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 2),
           tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 3));

    printf("\n");

    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 0),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 1),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 2),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 0),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 1),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 2),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 0),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 1),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 2),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e",
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 0),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 1),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 2),
           tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 3));

    printf("\n");
}

*/

// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_parallel(const double y[], double complex f_u[], double fvector[],
                double complex f_u_vector[]) {
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];
    int i, j, k; // Einstein summation over indices v and w

    LOOP_ijk gamma_udd[i][j][k] = 0.;

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double complex A_u[4] = {0., 0., 0., 0.}; // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
    connection_udd(X_u, gamma_udd);
    // connection_num_udd(X_u, gamma_udd);

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

double complex inner_product_complex_complex(const double *X_u,
                                             double complex *A_u,
                                             double complex *B_u) {
    // Obtain the covariant metric at X_u
    double g_dd[4][4];
    int i, j;
    LOOP_ij g_dd[i][j] = 0.;
    metric_dd(X_u, g_dd);

    // Compute the dot produt
    double complex dotproduct = 0.;
    LOOP_ij dotproduct += g_dd[i][j] * A_u[i] * B_u[j];

    return dotproduct;
}

double complex inner_product_real_complex(const double *X_u, double *A_u,
                                          double complex *B_u) {
    // Obtain the covariant metric at X_u
    double g_dd[4][4];
    int i, j;
    LOOP_ij g_dd[i][j] = 0.;
    metric_dd(X_u, g_dd);

    // Compute the dot produt
    double complex dotproduct = 0.;
    LOOP_ij dotproduct += g_dd[i][j] * A_u[i] * B_u[j];

    return dotproduct;
}

// TODO WARNING: handle the case when p=0
void stokes_to_f_tetrad_u(const double complex S_A[4],
                          double complex f_tetrad_u[4], double *p) {

    *p = 1.;
    //    if(1 || cabs(S_A[0]) > 1.e-150)
    //       *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) /
    //       S_A[0];

    if (cabs(S_A[1] + S_A[2] + S_A[3]) <
        1.e-100) { // Taylor expand for small values
        *p = (S_A[1] + S_A[2] + S_A[3]) / S_A[0];
    } else {
        *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) / S_A[0];
    }

    double complex QJones = S_A[1] / (S_A[0] * *p);
    double complex UJones = S_A[2] / (S_A[0] * *p);
    double complex VJones = S_A[3] / (S_A[0] * *p);

    // source:
    // https://physics.stackexchange.com/questions/238957/converting-stokes-parameters-to-jones-vector
    f_tetrad_u[1] = sqrt((1. + QJones) / 2.);

    if (f_tetrad_u[1] == 0)
        f_tetrad_u[2] = 1.;
    else
        f_tetrad_u[2] =
            UJones / (2. * f_tetrad_u[1]) - I * VJones / (2. * f_tetrad_u[1]);

    f_tetrad_u[1] *= sqrt(S_A[0] * *p);
    f_tetrad_u[2] *= sqrt(S_A[0] * *p);

    // A problem can occur where p is zero.
    // If p is zero, we have no polarization; vector is then irrelevant.
    if (*p < 1.e-30 || f_tetrad_u[1] != f_tetrad_u[1] ||
        f_tetrad_u[2] != f_tetrad_u[2]) {
        f_tetrad_u[1] = 0.;
        f_tetrad_u[2] = 1.;
        *p = 1.e-30;
    }

    if (f_tetrad_u[1] != f_tetrad_u[1] || f_tetrad_u[2] != f_tetrad_u[2]) {
        printf("PACKING NANS INTO VECTOR.\n");
        printf("\np = %+.15e", *p);
    }
}

void f_to_stokes(double Iinv, double Iinv_pol, double complex f_tetrad_u[],
                 double complex S_A[]) {
    S_A[0] = Iinv;
    S_A[1] = Iinv_pol * (cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) -
                         cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]));
    S_A[2] = Iinv_pol * (conj(f_tetrad_u[1]) * f_tetrad_u[2] +
                         f_tetrad_u[1] * conj(f_tetrad_u[2]));
    S_A[3] = Iinv_pol * (I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] -
                              f_tetrad_u[1] * conj(f_tetrad_u[2])));
}

void stokes_to_f(double complex S_A[], double *Iinv, double *Iinv_pol,
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

void f_tetrad_u_to_stokes(const double complex f_tetrad_u[4], const double p,
                          double complex S_A[4]) {
    double complex IfromJones = 0.;
    double complex QfromJones =
        0.; // cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) -
            // cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]);
    double complex UfromJones = 0.; // conj(f_tetrad_u[1]) * f_tetrad_u[2] +
                                    // f_tetrad_u[1] * conj(f_tetrad_u[2]);
    double complex VfromJones = 0.; // I * (conj(f_tetrad_u[1]) * f_tetrad_u[2]
                                    // - f_tetrad_u[1] * conj(f_tetrad_u[2]));

    IfromJones = (cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) +
                  cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2])) /
                 p;
    QfromJones = cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) -
                 cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]);
    UfromJones = conj(f_tetrad_u[1]) * f_tetrad_u[2] +
                 f_tetrad_u[1] * conj(f_tetrad_u[2]);
    VfromJones = I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] -
                      f_tetrad_u[1] * conj(f_tetrad_u[2]));

    S_A[0] = IfromJones;
    S_A[1] = QfromJones;
    S_A[2] = UfromJones;
    S_A[3] = VfromJones;

    if (S_A[0] != S_A[0] || S_A[1] != S_A[1] || S_A[2] != S_A[2] ||
        S_A[3] != S_A[3]) {
        printf("\n AND P IS %+.15e", p);
        printf("\nf_u[1] = %+.15e, %+.15e; f_u[2] = %+.15e, %+.15e",
               creal(f_tetrad_u[1]), cimag(f_tetrad_u[1]), creal(f_tetrad_u[2]),
               cimag(f_tetrad_u[2]));
        printf("NANS IN DA HOUSE \n NANS IN DA HOUSE \n NANS IN DA HOUSE");
    }
}

// NOTE: works only in Kerr metric
// Ziri's suggestion: construct U vecs
void construct_U_vector(const double X_u[], double U_u[]) {
    // Obtain relevant metric terms:
    double g_uu[4][4];
    metric_uu(X_u, g_uu);
    double g_uu00 = g_uu[0][0];
    double g_uu03 = g_uu[0][3];
    double g_uu33 = g_uu[3][3];

    // Observer/plasma wave vector:
    double U_d[4] = {-1., 0., 0., 0.};
    double B__ = -g_uu03 * U_d[0] / g_uu33;
    double C__ = -(1. + g_uu00 * U_d[0] * U_d[0]) / g_uu33;

    // Properly normalize U_u:
    U_d[3] = B__ + sqrt(B__ * B__ + C__);
    int i;
    LOOP_i U_u[i] = 0.;
    raise_index(X_u, U_d, U_u);

    // Check that U dot U = -1 to (near) machine precision:
    //    printf("\nU dot U: %+.15e", inner_product(X_u, U_u, U_u));
}

/*
#define ELECTRON_CHARGE    (4.80320425e-10)
#define ELECTRON_MASS      (9.1093829e-28)
#define PROTON_MASS        (1.6726219e-24)
#define BOLTZMANN_CONSTANT (1.3806488e-16)
#define SPEED_OF_LIGHT     (2.99792458e10)
#define PLANCK_CONSTANT    (6.62606885e-27)
#define MPCL2              (0.0015033)
#define GGRAV              (6.674e-8)
#define MSUN               (1.989e33)
#define MPoME              (PROTON_MASS / ELECTRON_MASS)
#define M_PI           3.14159265358979323846
*/

// Dexter (2016) A.18
double I_I(double x) {
    return 2.5651 * (1. + 1.92 * pow(x, -1. / 3.) + 0.9977 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.19
double I_Q(double x) {
    return 2.5651 *
           (1. + 0.932 * pow(x, -1. / 3.) + 0.4998 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.20
double I_V(double x) {
    return (1.8138 / x + 3.423 * pow(x, -2. / 3.) + 0.02955 * pow(x, -0.5) +
            2.0377 * pow(x, -1 / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double I_I_(double x) {
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) + 0.9977 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double I_Q_(double x) {
    return 2.5651 *
           (1 + 0.93193 * pow(x, -1. / 3.) + 0.499873 * pow(x, -2. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

double I_V_(double x) {
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
            0.0292545 * pow(x, -0.5) + 2.03773 * pow(x, -1. / 3.)) *
           exp(-1.8899 * pow(x, 1. / 3.));
}

// A.12
double j_I(double theta_e, double n_e, double nu, double B, double theta_B) {

    // double nu_B = e * B / (2. * M_PI * m * c);
    // double nu_p = 3. / 2. * nu_B * sin(theta_B);
    // double nu_c = nu_p * theta_e * theta_e;
    // Alternatively,
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    //    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / (2. * sqrt(3.) *
    //    SPEED_OF_LIGHT * theta_e * theta_e) * I_I(x);

    // NEW VERSION
    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_I_(x);
}

// A.13
double j_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    // return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / (2. * sqrt(3.) *
    // SPEED_OF_LIGHT * theta_e * theta_e) * I_Q(x);

    // NEW VERSION
    return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           SPEED_OF_LIGHT / theta_e / theta_e * I_Q_(x);
}

// A.14 (theta_B = pitch angle, k dot B)
double j_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double x = nu / nu_c;

    //    return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu *
    //    tan(theta_B) / (3. * sqrt(3.) * SPEED_OF_LIGHT * theta_e * theta_e *
    //    theta_e) * I_V(x);

    // NEW VERSION
    return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / tan(theta_B) /
           3. / sqrt(3.) / SPEED_OF_LIGHT / theta_e / theta_e / theta_e *
           I_V_(x);
}

// ABSORPTION
/////////////

// Use Kirchoff: a_nu = j_nu/B_nu

// ROTATION
///////////

// B.6
double f(double X) {
    return 2.011 * exp(-pow(X, 1.035) / 4.7) -
           cos(X / 2.) * exp(-sqrt(X) / 2.73) - 0.011 * exp(-X / 47.2);
}

// B.13 NOTE: make sure that, in B.13, Dexter does not mean X when he writes x
// (in the final term) Note that X = ... (B.8) UPDATE: it was apparently indeed
// a typo. while x = nu / nu_c UPDATE: it seems Dexter did make TWO typos: it
// should be ln(X/120), not "ln x / 120". Thus I have greyed out the original,
// new one is below.
double f_m(double X, double x) {
    //    return f(X) + (0.011 * exp(-X/47.2) - pow(2.,-1./3.)/pow(3.,23./6.) *
    //    10000. * M_PI * pow(X, -8./3.)) * 0.5 * (1. + tanh(10. *
    //    log(x)/120.));
    return f(X) +
           (0.011 * exp(-X / 47.2) - pow(2., -1. / 3.) / pow(3., 23. / 6.) *
                                         10000. * M_PI * pow(X, -8. / 3.)) *
               0.5 * (1. + tanh(10. * log(X / 120.)));
}

// Bessel function approximations:
double K_0(double x) { return -log(0.5 * x) - 0.5772; }
double K_1(double x) { return 1. / x; }
double K_2(double x) { return 2. / (x * x); }

// B.15:
double DeltaJ_5(double X) {
    return 0.4379 * log(1. + 0.001858 * pow(X, 1.503));
}

/*
    *rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu,
   4) * jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
         6. * Thetae) * sin(theta) * sin(theta);
    *rU = 0.0;
    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
        (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) *
   cos(theta);
*/

double besselk_asym(int n, double x) {

    //    return gsl_sf_bessel_Kn(n, x);

    if (n == 0)
        return -log(x / 2.) - 0.5772;

    if (n == 1)
        return 1. / x;

    if (n == 2)
        return 2. / x / x;

    fprintf(stderr, "this cannot happen\n");
    exit(1);
}

double jffunc(double Xe) {
    double extraterm;
    extraterm =
        (0.011 * exp(-Xe / 47.2) - pow(2., -1. / 3.) / pow(3., 23. / 6.) *
                                       M_PI * 1e4 * pow(Xe + 1e-16, -8. / 3.)) *
        (0.5 + 0.5 * tanh((log(Xe) - log(120.)) / 0.1));
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
           cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) - 0.011 * exp(-Xe / 47.2) +
           extraterm;
}

double Je(double Xe) {
    return 0.43793091 * log(1. + 0.00185777 * pow(Xe, 1.50316886));
}

// B.4
double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_B =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double X = sqrt(3. / (2. * sqrt(2.)) * 0.001 * nu / nu_c);
    double x = nu / nu_c;

    //  return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu_B * nu_B *
    //  sin(theta_B) * sin(theta_B) / (ELECTRON_MASS * SPEED_OF_LIGHT * nu * nu
    //  * nu) *
    //       f_m(X,x) * (K_1(1./theta_e) / K_2(1./theta_e) + 6. * theta_e);

    // NEW VERSION
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    return 2. * M_PI * nu / 2. / SPEED_OF_LIGHT * wp2 * omega0 * omega0 /
           pow(2. * M_PI * nu, 4.) * jffunc(Xe) *
           (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
            6. * theta_e) *
           sin(theta_B) * sin(theta_B);
}

// B.14
double rho_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    double nu_B =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e +
                  1.0;

    double X = sqrt(3. / (2. * sqrt(2.)) * 0.001 * nu / nu_c);

    //  return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu_B /
    //  (ELECTRON_MASS * SPEED_OF_LIGHT * nu * nu) * K_0(1. / theta_e) -
    //  DeltaJ_5(X) /
    //        K_2(1. / theta_e);

    // NEW VERSION
    double wp2 =
        4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    double Thetaer = 1. / theta_e;

    return 2.0 * M_PI * nu / SPEED_OF_LIGHT * wp2 * omega0 /
           pow(2. * M_PI * nu, 3) * (besselk_asym(0, Thetaer) - Je(Xe)) /
           besselk_asym(2, Thetaer) * cos(theta_B);
}

double radiative_transfer(double *lightpath, int steps, double frequency) {
    int IN_VOLUME, path_counter;
    double I_current = 0.;
    double dI = 0.;
    double j_nu = 0.;
    double B, THETA_e, pitch_ang, nu_p, n_e, nu_p2, dl_current;
    int i;
    double X_u[4], k_u[4], k_d[4], B_u[4], Uplasma_u[4];
    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    double tau = 0.;
    double a_nu = 0.;
    double K_inv_old = 0, j_inv_old = 0, dtau_old = 0;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        // Obtain the parameters n_e, THETA_e, B, and Uplasma_u at X_u
        // get_plasma_parameters(X_u, &n_e, &THETA_e, &B, Uplasma_u);
        get_fluid_params(X_u, &n_e, &THETA_e, &B, B_u, Uplasma_u, &IN_VOLUME);

        // Check whether the ray is currently in the GRMHD simulation volume
        if (IN_VOLUME) {
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

            // lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
            nu_p2 = nu_p * nu_p;

            // Obtain emission coefficient in current plasma conditions
            j_nu = emission_coeff_THSYNCHAV(B, THETA_e, nu_p, n_e);

            // Obtain absorption coefficient
            if (ABSORPTION) {
                a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);
            }

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT /
                       (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            double redshift = frequency / nu_p;

            double dtau = (nu_p * a_nu * dl_current * C + dtau_old);
            double K_inv = (nu_p * a_nu);
            double j_inv = (j_nu / nu_p2);

            // Only add I_current if it is not NaN
            if (j_nu == j_nu &&
                exp(X_u[1]) <
                    RT_OUTER_CUTOFF) { // I_current += exp(-tau) * j_nu /
                                       // nu_p / nu_p * dl_current * C;
                //         I_current += dI; // Old way of integrating
                double Ii = I_current;
                double S = j_inv / K_inv;
                if (K_inv == 0)
                    I_current = Ii;
                else if (dtau < 1.e-5)
                    I_current = Ii - (Ii - S) * (0.166666667 * dtau *
                                                 (6. - dtau * (3. - dtau)));
                else {
                    double efac = exp(-dtau);
                    I_current = Ii * efac + S * (1. - efac);
                }
            }
        }
    }

    // Store integrated intensity in the image
    return I_current * pow(frequency, 3.);
}

// TODO WARNING: handle the case when p=0
void stokes_to_f_tetrad_u_old(const double complex S_A[4],
                              double complex f_tetrad_u[4], double *p) {
    // *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) / S_A[0];

    if (cabs(S_A[1] + S_A[2] + S_A[3]) <
        1.e-100) { // Taylor expand for small values
        *p = (S_A[1] + S_A[2] + S_A[3]) / S_A[0];
    } else {
        *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) / S_A[0];
    }

    double complex QJones = S_A[1] / (S_A[0] * *p);
    double complex UJones = S_A[2] / (S_A[0] * *p);
    double complex VJones = S_A[3] / (S_A[0] * *p);

    // source:
    // https://physics.stackexchange.com/questions/238957/converting-stokes-parameters-to-jones-vector

    f_tetrad_u[1] = sqrt((1. + QJones) / 2.);

    if (f_tetrad_u[1] == 0)
        f_tetrad_u[2] = 1.;
    else
        f_tetrad_u[2] =
            UJones / (2. * f_tetrad_u[1]) - I * VJones / (2. * f_tetrad_u[1]);

    f_tetrad_u[1] *= sqrt(S_A[0] * *p);
    f_tetrad_u[2] *= sqrt(S_A[0] * *p);

    //    fprintf(stderr, "\n\nf dot f %+.15e\n\n", f_tetrad_u[1] *
    //    conj(f_tetrad_u[1]) + f_tetrad_u[2] * conj(f_tetrad_u[2]));
}

double check_handedness(double X_u[4], double tetrad_u[][4]) {
    int i, j, k, l;
    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double g = determ(g_dd, 4);

    double eps[4][4][4][4];
    LOOP_ijkl {
        if ((i == j) || (i == k) || (i == l) || (j == k) || (j == l) ||
            (k == l))
            eps[i][j][k][l] = 0.;
        else
            eps[i][j][k][l] = ((i - j) * (i - k) * (i - l) * (j - k) * (j - l) *
                               (k - l) / 12.);
    }

    double result = 0.;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            for (l = 0; l < 4; l++)
                for (k = 0; k < 4; k++) {
                    result += g * eps[i][j][k][l] * tetrad_u[i][0] *
                              tetrad_u[j][1] * tetrad_u[k][2] * tetrad_u[l][3];
                }

    return result;
}

double radiative_transfer_polarized(double *lightpath, int steps,
                                    double frequency, double *f_x, double *f_y,
                                    double *p, int PRINT_POLAR, double *IQUV) {
    int IN_VOLUME, path_counter;
    double I_current = 0.;
    double dI = 0.;
    double j_nu = 0.;
    double B, THETA_e, pitch_ang, nu_p, n_e, nu_p2, dl_current;
    int i, j;
    double X_u[4], k_u[4], k_d[4], B_u[4], Uplasma_u[4];
    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    double tau = 0.;
    double a_nu = 0.;
    double K_inv_old = 0, j_inv_old = 0, dtau_old = 0;

    int POLARIZATION_ACTIVE = 0;
    double S_I_current = 0.;
    double S_Q_current = 0.;
    double S_U_current = 0.;
    double S_V_current = 0.;
    double S_I_new, S_Q_new, S_U_new, S_V_new;
    double deg_of_pol = 0.;

    double tetrad_u[4][4], tetrad_d[4][4];
    LOOP_ij tetrad_u[i][j] = 0.;
    LOOP_ij tetrad_d[i][j] = 0.;

    double photon_u_current[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

    double complex f_tetrad_u[4] = {0., 0., 0., 0.};
    double complex f_u[4] = {0., 0., 0., 0.};

    double complex S_A[4] = {0., 0., 0., 0.};

    double Iinv = 0.;
    double Iinv_pol = 0.;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        get_fluid_params(X_u, &n_e, &THETA_e, &B, B_u, Uplasma_u, &IN_VOLUME);

        // PLASMA INTEGRATION STEP
        //////////////////////////

        double r_current2 = logscale ? exp(X_u[1]) : X_u[1];
        double OUTER_BOUND_POL = 1000.;

        // Check whether the ray is currently in the GRMHD simulation volume
        if (IN_VOLUME && r_current2 < OUTER_BOUND_POL) {
            // Obtain pitch angle: still no units (geometric)
            pitch_ang = pitch_angle(X_u, k_u, B_u, Uplasma_u);

            // CGS UNITS USED FROM HERE ON OUT
            //////////////////////////////////

            // For the polarized case, we have to do it differently.
            // Unpolarized: 1) Create light path by integration. 2) For each
            // step in lightpath, perform one radiative transfer step.
            // Polarized:   1) Create light path by integration. 2) For each
            // step in lightpath, perform one radiative transfer step, AND,
            // OUTSIDE in_volume loop, do a spacetime propagation step.

            // TRANSFER STEP
            ////////////////

            // Scale the wave vector to correct energy
            LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                             (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // Convert distance dlambda accordingly
            dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                          (PLANCK_CONSTANT * frequency);

            // lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
            nu_p2 = nu_p * nu_p;

            // UNPOLARIZED EMISSION/ABSORPTION COEFFS
            /////////////////////////////////////////

            // Obtain emission coefficient in current plasma conditions
            j_nu = emission_coeff_THSYNCHAV(B, THETA_e, nu_p, n_e);

            // Obtain absorption coefficient
            if (ABSORPTION) {
                a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);
            }

            // POLARIZED EMISSION/ABSORPTION COEFFS
            ///////////////////////////////////////

            double jI = j_I(THETA_e, n_e, nu_p, B, pitch_ang);
            double jQ = j_Q(THETA_e, n_e, nu_p, B, pitch_ang);
            double jU = 0.;
            double jV = j_V(THETA_e, n_e, nu_p, B, pitch_ang);

            double rQ = rho_Q(THETA_e, n_e, nu_p, B, pitch_ang);
            double rU = 0.;
            double rV = rho_V(THETA_e, n_e, nu_p, B, pitch_ang);

            double aI = absorption_coeff_TH(jI, nu_p, THETA_e);
            double aQ = absorption_coeff_TH(jQ, nu_p, THETA_e);
            double aU = absorption_coeff_TH(jU, nu_p, THETA_e);
            double aV = absorption_coeff_TH(jV, nu_p, THETA_e);

            // Transform to invariant forms
            jI /= (nu_p * nu_p);
            jQ /= (nu_p * nu_p);
            jU /= (nu_p * nu_p);
            jV /= (nu_p * nu_p);

            aI *= nu_p;
            aQ *= nu_p;
            aU *= nu_p;
            aV *= nu_p;

            rQ *= nu_p;
            rU *= nu_p;
            rV *= nu_p;

            // UNPOLARIZED TRANSFER STEP
            ////////////////////////////

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT /
                       (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            double redshift = frequency / nu_p;

            j_nu = j_I(THETA_e, n_e, nu_p, B, pitch_ang);
            a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);

            double K_inv = aI;
            double dtau = (K_inv * dl_current * C + dtau_old);
            double j_inv = jI;

            // POLARIZED TRANSFER
            /////////////////////

            // Create tetrad, needed whether POLARIZATION_ACTIVE is true or
            // false.
            create_observer_tetrad_u2(X_u, k_u, Uplasma_u, B_u, tetrad_u);
            create_tetrad_d(X_u, tetrad_u, tetrad_d);

            double tetradCheck = check_tetrad_compact(X_u, tetrad_u);

            if (tetradCheck != tetradCheck)
                printf("NAN TETRAD\n");

            // FROM F VECTOR TO STOKES (when applicable)
            ////////////////////////////////////////////

            // If (POLARIZATION_ACTIVE), get Stokes params from f_u and p.
            // (Otherwise, never been in volume before; we simply use
            // S_I_current)
            if (POLARIZATION_ACTIVE && 0) {
                // Get f_tetrad_u
                LOOP_i f_tetrad_u[i] = 0.;
                LOOP_ij f_tetrad_u[i] += tetrad_d[j][i] * f_u[j];

                // Get Stokes params from f_tetrad_u
                f_tetrad_u_to_stokes(f_tetrad_u, *p, S_A);
            }

            // NEW VERSION
            if (POLARIZATION_ACTIVE && 1) {
                // Get f_tetrad_u
                LOOP_i f_tetrad_u[i] = 0.;
                LOOP_ij f_tetrad_u[i] += tetrad_d[j][i] * f_u[j];

                // Get Stokes params from f_tetrad_u
                f_to_stokes(Iinv, Iinv_pol, f_tetrad_u, S_A);
            }

            // Given Stokes params and plasma coeffs, compute NEW Stokes params
            // after plasma step.
            double complex I0 = S_A[0];
            double complex Q0 = S_A[1];
            double complex U0 = S_A[2];
            double complex V0 = S_A[3];

            double THRESH = 0.1;

            // New stiffness check
            double a2 = rQ * rQ + rV * rV - aQ * aQ - aV * aV;
            double a0 =
                -2. * aV * aQ * rV * rQ - aQ * aQ * rQ * rQ - aV * aV * rV * rV;

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

            if (res1 > STIFFTHRESH || res2 > STIFFTHRESH ||
                res3 > STIFFTHRESH || res4 > STIFFTHRESH)
                STIFF = 1;

            // If both rotation coeffs (times dlambda) are smaller than
            // threshold, take an RK4 step; otherwise, implicit Euler.
            // if(!STIFF)
            if (fabs(rQ) < THRESH && fabs(rV) < THRESH) {
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
                double complex Ik2 =
                    dl_current * C * jI -
                    dl_current * C *
                        (aI * (I0 + 0.5 * Ik1) + aQ * (Q0 + 0.5 * Qk1) +
                         aU * (U0 + 0.5 * Uk1) + aV * (V0 + 0.5 * Vk1));
                double complex Qk2 =
                    dl_current * C * jQ -
                    dl_current * C *
                        (aQ * (I0 + 0.5 * Ik1) + aI * (Q0 + 0.5 * Qk1) +
                         rV * (U0 + 0.5 * Uk1) - rU * (V0 + 0.5 * Vk1));
                double complex Uk2 =
                    dl_current * C * jU -
                    dl_current * C *
                        (aU * (I0 + 0.5 * Ik1) - rV * (Q0 + 0.5 * Qk1) +
                         aI * (U0 + 0.5 * Uk1) + rQ * (V0 + 0.5 * Vk1));
                double complex Vk2 =
                    dl_current * C * jV -
                    dl_current * C *
                        (aV * (I0 + 0.5 * Ik1) + rU * (Q0 + 0.5 * Qk1) -
                         rQ * (U0 + 0.5 * Uk1) + aI * (V0 + 0.5 * Vk1));

                // k3
                double complex Ik3 =
                    dl_current * C * jI -
                    dl_current * C *
                        (aI * (I0 + 0.5 * Ik2) + aQ * (Q0 + 0.5 * Qk2) +
                         aU * (U0 + 0.5 * Uk2) + aV * (V0 + 0.5 * Vk2));
                double complex Qk3 =
                    dl_current * C * jQ -
                    dl_current * C *
                        (aQ * (I0 + 0.5 * Ik2) + aI * (Q0 + 0.5 * Qk2) +
                         rV * (U0 + 0.5 * Uk2) - rU * (V0 + 0.5 * Vk2));
                double complex Uk3 =
                    dl_current * C * jU -
                    dl_current * C *
                        (aU * (I0 + 0.5 * Ik2) - rV * (Q0 + 0.5 * Qk2) +
                         aI * (U0 + 0.5 * Uk2) + rQ * (V0 + 0.5 * Vk2));
                double complex Vk3 =
                    dl_current * C * jV -
                    dl_current * C *
                        (aV * (I0 + 0.5 * Ik2) + rU * (Q0 + 0.5 * Qk2) -
                         rQ * (U0 + 0.5 * Uk2) + aI * (V0 + 0.5 * Vk2));

                // k4
                double complex Ik4 = dl_current * C * jI -
                                     dl_current * C *
                                         (aI * (I0 + Ik3) + aQ * (Q0 + Qk3) +
                                          aU * (U0 + Uk3) + aV * (V0 + Vk3));
                double complex Qk4 = dl_current * C * jQ -
                                     dl_current * C *
                                         (aQ * (I0 + Ik3) + aI * (Q0 + Qk3) +
                                          rV * (U0 + Uk3) - rU * (V0 + Vk3));
                double complex Uk4 = dl_current * C * jU -
                                     dl_current * C *
                                         (aU * (I0 + Ik3) - rV * (Q0 + Qk3) +
                                          aI * (U0 + Uk3) + rQ * (V0 + Vk3));
                double complex Vk4 = dl_current * C * jV -
                                     dl_current * C *
                                         (aV * (I0 + Ik3) + rU * (Q0 + Qk3) -
                                          rQ * (U0 + Uk3) + aI * (V0 + Vk3));

                S_A[0] = I0 + 1. / 6. * (Ik1 + 2. * Ik2 + 2. * Ik3 + Ik4);
                S_A[1] = Q0 + 1. / 6. * (Qk1 + 2. * Qk2 + 2. * Qk3 + Qk4);
                S_A[2] = U0 + 1. / 6. * (Uk1 + 2. * Uk2 + 2. * Uk3 + Uk4);
                S_A[3] = V0 + 1. / 6. * (Vk1 + 2. * Vk2 + 2. * Vk3 + Vk4);
            } else {
                /*
                                double u11 = 1. + dl_current * C * aI;
                                double u12 = dl_current * C * aQ;
                                double u14 = dl_current * C * aV;
                                double l21 = dl_current * C * aQ / u11;
                                double u22 = 1. + dl_current * C * aI - l21 *
                   u12; double u23 = dl_current * C * rV; double u24 = -l21 *
                   u14; double l32 = -dl_current * C * rV / u22; double u33 = 1.
                   + dl_current * C * aI - l32 * u23; double u34 = dl_current *
                   C * rQ - l32 * u24; double l41 = dl_current * C * aV / u11;
                                double l42 = -l41 * u12 / u22;
                                double l43 = (-dl_current * C * rQ - l42 * u23)
                   / u33; double u44 = 1. + dl_current * C * aI - l41 * u14 -
                   l42
                   * u24 - l43 * u34;
                */
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
                double u44 = 1. + 0.5 * dl_current * C * aI - l41 * u14 -
                             l42 * u24 - l43 * u34;

                // Construct b-vector.
                /*
                                double b1 = I0 + dl_current * C * jI;
                                double b2 = Q0 + dl_current * C * jQ;
                                double b3 = U0 + dl_current * C * jU;
                                double b4 = V0 + dl_current * C * jV;
                */
                double b1 = I0 + dl_current * C / 2. *
                                     (2. * jI - (aI * I0 + aQ * Q0 + aV * V0));
                double b2 = Q0 + dl_current * C / 2. *
                                     (2. * jQ - (aQ * I0 + aI * Q0 + rV * U0));
                double b3 = U0 + dl_current * C / 2. *
                                     (2. * jU - (-rV * Q0 + aI * U0 + rQ * V0));
                double b4 = V0 + dl_current * C / 2. *
                                     (2. * jV - (aV * I0 - rQ * U0 + aI * V0));

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

            // FROM STOKES TO F VECTOR
            ///////////////////////////

            Iinv = S_A[0];

            Iinv_pol =
                sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

            double Iinvt = 0.;
            double Iinv_polt = 0.;
            double S_Atest[4] = {0., 0., 0., 0.};
            double f_test[4] = {0., 0., 0., 0.};
            if (Iinv_pol > 1.e-70 && 0) {
                fprintf(stderr, "\n\nPRE: %g, %g, %g, %g ", S_A[0], S_A[1],
                        S_A[2], S_A[3]);
                stokes_to_f(S_A, &Iinvt, &Iinv_polt, f_test);
                fprintf(stderr, " Iinvt = %g ", Iinvt);
                f_to_stokes(Iinvt, Iinv_polt, f_test, S_Atest);
                fprintf(stderr, "POST: %g, %g, %g, %g \n\n", S_Atest[0],
                        S_Atest[1], S_Atest[2], S_Atest[3]);
            }

            if (Iinv_pol > 1.e-100 && 1) {
                stokes_to_f(S_A, &Iinv, &Iinv_pol, f_tetrad_u);

                // Update f_u using f_tetrad_u.
                LOOP_i f_u[i] = 0.;
                LOOP_ij f_u[i] += tetrad_u[i][j] * f_tetrad_u[j];

                // Set POLARIZATION_ACTIVE to true; we are, after all,
                // in_volume.
                POLARIZATION_ACTIVE = 1;

            } else {
                POLARIZATION_ACTIVE = 0;
                S_A[1] = 0.;
                S_A[2] = 0.;
                S_A[3] = 0.;
            }

            // We have now updated the Stokes vector using plasma at current
            // position. Only do stuff below this line IF S_A[0] > 1.e-40. If
            // not, POLARIZATION_ACTIVE is set to FALSE and we reset S_A[i] = 0
            double STOKES_I_THRESHOLD = 1.e-70;
            if (cabs(S_A[0]) > STOKES_I_THRESHOLD && 0) {
                if (cabs(S_A[1] + S_A[2] + S_A[3]) <
                    1.e-100) { // Taylor expand for small values
                    *p = (S_A[1] + S_A[2] + S_A[3]) / S_A[0];
                } else {
                    *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] +
                              S_A[3] * S_A[3]) /
                         S_A[0];
                }

                // Pack new Stokes params into f_u and p.
                stokes_to_f_tetrad_u(S_A, f_tetrad_u, p);

                // Update f_u using f_tetrad_u.
                LOOP_i f_u[i] = 0.;
                LOOP_ij f_u[i] += tetrad_u[i][j] * f_tetrad_u[j];

                // Set POLARIZATION_ACTIVE to true; we are, after all,
                // in_volume.
                POLARIZATION_ACTIVE = 1;

                // Get f_tetrad_u
                // LOOP_i  f_tetrad_u[i] = 0.;
                // LOOP_ij f_tetrad_u[i] += tetrad_d[j][i] * f_u[j];

            } // else{
              //    POLARIZATION_ACTIVE = 0;
              //    LOOP_i S_A[i] = 0.;
              //}
        }     // End of if(IN_VOLUME)

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

            // Obtain the right dlambda.
            // Already good.

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
    double cam_up_u[4] = {0., 0., 0., -1.};

    if (0) {
        printf("\n X_u[0] = %+.15e", X_u[0]);
        printf("\n X_u[1] = %+.15e", X_u[1]);
        printf("\n X_u[2] = %+.15e", X_u[2]);
        printf("\n X_u[3] = %+.15e", X_u[3]);
    }

    // Need U_obs_u
    double U_obs_u[4] = {0., 0., 0., 0.};
    double obs_tetrad_u[4][4], obs_tetrad_d[4][4];
    LOOP_ij obs_tetrad_u[i][j] = 0.;
    LOOP_ij obs_tetrad_d[i][j] = 0.;
    construct_U_vector(X_u, U_obs_u);

    create_observer_tetrad_u2(X_u, k_u, U_obs_u, cam_up_u, obs_tetrad_u);
    create_tetrad_d(X_u, obs_tetrad_u, obs_tetrad_d);

    double tetradCheck = check_tetrad_compact(X_u, obs_tetrad_u);

    // printf("observer handedness: %g\n", check_handedness(X_u, obs_tetrad_u));

    //    printf("\n Observer tetradcheck %+.15e", tetradCheck);
    //    printf("\n Observer U dot U %+.15e", inner_product(X_u, U_obs_u,
    //    U_obs_u)); printf("\n Observer k dot k %+.15e", inner_product(X_u,
    //    k_u, k_u)); printf("\n Observer k_rev dot k_rev %+.15e",
    //    inner_product(X_u, k_u_reverse, k_u_reverse));

    // Convert f_u to f_obs_tetrad_u
    double complex f_obs_tetrad_u[4] = {0., 0., 0., 0.};
    LOOP_i f_obs_tetrad_u[i] = 0.;
    LOOP_ij f_obs_tetrad_u[i] += obs_tetrad_d[j][i] * f_u[j];

    double complex S_If = 0.;
    double complex S_Qf = 0.;
    double complex S_Uf = 0.;
    double complex S_Vf = 0.;

    if (POLARIZATION_ACTIVE && 0) {
        f_tetrad_u_to_stokes(f_obs_tetrad_u, *p, S_A);

        // Construct final (NON-INVARIANT) Stokes params.
        S_If = S_A[0] * pow(frequency, 3.);
        S_Qf = S_A[1] * pow(frequency, 3.);
        S_Uf = S_A[2] * pow(frequency, 3.);
        S_Vf = S_A[3] * pow(frequency, 3.);
    }

    if (POLARIZATION_ACTIVE && 1) {
        f_to_stokes(Iinv, Iinv_pol, f_obs_tetrad_u, S_A);

        // Construct final (NON-INVARIANT) Stokes params.
        S_If = S_A[0] * pow(frequency, 3.);
        S_Qf = S_A[1] * pow(frequency, 3.);
        S_Uf = S_A[2] * pow(frequency, 3.);
        S_Vf = S_A[3] * pow(frequency, 3.);
    }

    IQUV[0] = S_If;
    IQUV[1] = S_Qf;
    IQUV[2] = S_Uf;
    IQUV[3] = S_Vf;

    if (cabs(S_If) > 1.e-40 && 0) {
        printf("\nS_If = %+.15e", creal(S_If));
        printf("\nS_Qf = %+.15e", creal(S_Qf));
        printf("\nS_Uf = %+.15e", creal(S_Uf));
        printf("\nS_Vf = %+.15e", creal(S_Vf));
        *p = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) / S_A[0];
        printf("\npf = %+.15e", *p);
    }

    // Read final polarization variables at observer.
    *f_x = 0.;
    *f_y = 0.;
    *p = 0.;
    //    *deg_of_pol = 0.;

    // Stokes...

    // Store integrated intensity in the image.
    return S_If; // I_current * pow(frequency, 3.);
}
