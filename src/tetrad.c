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

// FUNCTIONS
////////////

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

// Creates a tetrad with a random orientation, i.e., z-axis fixed but XY
// axes can have any azimuth.
// (Useful when we have no orienting vector and just want to express
// our polarized state in terms of Stokes params.)
void create_tetrad(double X_u[], double k_u[], double U_u[],
                   double tetrad_u[][4]) {
    // Summation indices:

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double e_u_t[4];
    LOOP_i e_u_t[i] = U_u[i];

    double e_u_K[4];
    double omega = -inner_product(X_u, k_u, U_u);
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i];

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
    LOOP_i b_u[i] = 1. / sqrt(3.) * (double)rand() / (double)RAND_MAX;

    // 2. Check inner product.
    // fprintf(stderr, "\nb dot b = %+.15e", inner_product(X_u, b_u, b_u));

    // 3. If spacelike, OK.
    // 4. If not, flip sign of b_u[0]? Or just create new b_u?
    while (inner_product(X_u, b_u, b_u) < 0.1) {
        LOOP_i b_u[i] = 1. / sqrt(3.) * (double)rand() / (double)RAND_MAX;
    }

    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);

    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] =
        (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive;

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

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] +=
        (-1. / sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) /
        (omega * Ncursive);

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

// Creates a tetrad whose Y-axis is aligned with vector b_u
// (useful when we have such an orienting vector, e.g., observer)
void create_observer_tetrad(double X_u[], double k_u[], double U_u[],
                            double Bs_u[], double tetrad_u[][4]) {
    // Summation indices:
    double b_u[4];
    LOOP_i b_u[i] = Bs_u[i] * B_unit;

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double e_u_t[4];
    LOOP_i e_u_t[i] = U_u[i];

    double e_u_K[4];
    double omega = -inner_product(X_u, k_u, U_u);
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i];

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

    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);

    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] =
        (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive;

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

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] +=
        (-1. / sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) /
        (omega * Ncursive);

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

double tetrad_identity_eta(double X_u[4], double tetrad_u[4][4], int a, int b) {
    double result = 0.;

    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    LOOP_ij result += g_dd[i][j] * tetrad_u[i][a] * tetrad_u[j][b];

    return result;
}

double tetrad_identity_g(double tetrad_u[][4], int mu, int nu) {

    double eta_uu[4][4] = {{-1., 0., 0., 0.},
                           {0., 1., 0., 0.},
                           {0., 0., 1., 0.},
                           {0., 0., 0., 1.}};

    double result = 0.;

    LOOP_ij result += eta_uu[i][j] * tetrad_u[mu][i] * tetrad_u[nu][j];

    return result;
}

double tetrad_identity_sum_latin(double tetrad_u[4][4], double tetrad_d[4][4],
                                 int mu, int nu) {

    double result = 0.;

    LOOP_i result += tetrad_u[mu][i] * tetrad_d[nu][i];

    return result;
}

double tetrad_identity_sum_greek(double tetrad_u[4][4], double tetrad_d[4][4],
                                 int a, int b) {

    double result = 0.;

    LOOP_i result += tetrad_u[i][a] * tetrad_d[i][b];

    return result;
}

void create_tetrad_d(double X_u[], double tetrad_u[][4], double tetrad_d[][4]) {
    double eta_minkowski[4][4] = {
        {-1., 0., 0., 0.},
        {0., 1., 0., 0.},
        {0., 0., 1., 0.},
        {0., 0., 0., 1.},
    };

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

double check_tetrad_compact(double X_u[], double tetrad_u[][4]) {
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
    fprintf(stderr, "\nTetrad identities (should be close to zero): %+.15e",
            result);
    return result;
}

void check_tetrad_identities(double X_u[], double tetrad_u[][4]) {
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
