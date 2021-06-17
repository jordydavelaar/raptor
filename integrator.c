/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 */

#include <math.h>
#include "functions.h"
#include "parameters.h"
#include "constants.h"
#include <complex.h>

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk4_step(double *y, void (*f)(double*, double*), double dt){
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
    for (q = 0; q < 4; q++){
        f(yshift, fvector); // Apply function f to current y to obtain fvector
        for (i = 0; i < DIM * 2; i++){
            dx[q * DIM * 2 + i] = dt * fvector[i]; // Use fvector to fill dx
            yshift[i] = y[i] + dx[q * DIM * 2 + i] * weights[q]; // Update y
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < DIM * 2; i++)
        y[i] = y[i] + 1. / 6. * (dx[0 * DIM * 2 + i] +
                                 dx[1 * DIM * 2 + i] * 2. +
                                 dx[2 * DIM * 2 + i] * 2. +
                                 dx[3 * DIM * 2 + i]);
}

// Updates the vector y (containing position/velocity) by one RK2 step.
void rk2_step(double *y, void (*f)(double*, double*), double dt){
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
    for (q = 0; q < 2; q++){
        f(yshift, fvector); // Apply function f to current y to obtain fvector
        for (i = 0; i < DIM * 2; i++){
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
void verlet_step(double *y, void (*f)(double*, double*), double dl){
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
        yshift[i + DIM] = yshift[i + DIM] + fvector[i+ DIM] * dl;
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
double stepsize(double X_u[4], double U_u[4]){
    double SMALL = 1.e-40;

    double dlx1  = STEPSIZE / (fabs(U_u[1]) + SMALL*SMALL);
    double dlx2  = STEPSIZE * fmin( X_u[2], M_PI - X_u[2])/ (fabs(U_u[2]) + SMALL*SMALL);
    double dlx3  = STEPSIZE / (fabs(U_u[3]) + SMALL*SMALL);

    double idlx1 = 1. / (fabs(dlx1) + SMALL*SMALL) ;
    double idlx2 = 1. / (fabs(dlx2) + SMALL*SMALL) ;
    double idlx3 = 1. / (fabs(dlx3) + SMALL*SMALL) ;

    return -1. / (idlx1 + idlx2 + idlx3) ;
}

// The function to be used by the integrator for GR geodesic calculations.
// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_geodesic(double *y, double *fvector){
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double A_u[4] = {0.,   0.,   0.,   0.  }; // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
    connection_udd(X_u, gamma_udd);
    //connection_num_udd(X_u, gamma_udd);

    // Compute 4-acceleration using the geodesic equation
    int i, j, k; // Einstein summation over indices v and w
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];

    // Update fvector
    LOOP_i {
        fvector[i]       = U_u[i];
        fvector[i + DIM] = A_u[i];
    }
}















// Recursively calculates the determinant of a matrix.
// Source: http://ideone.com/fork/92JF0O
double determ(double matrix[][4],int n)
{
    int p, h, k, i, j;
    double det = 0.;
    double temp[4][4];
    if(n==1) {
        return matrix[0][0];
    } else if(n==2) {
        det=(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]);
        return det;
    } else {
        for(p=0;p<n;p++) {
            h = 0;
            k = 0;
            for(i=1;i<n;i++)
                for( j=0;j<n;j++) {
                    if(j==p)
                        continue;
                    temp[h][k] = matrix[i][j];
                    k++;
                    if(k==n-1) {
                        h++;
                        k = 0;
                    }
                }
            det=det+matrix[0][p]*pow(-1,p)*determ(temp,n-1);
        }
        return det;
    }
}


// Here we implement Ziri's tetrad
void create_tetrad_u2(const double X_u[], const double k_u[], const double U_u[], double tetrad_u[][4])
{
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
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i]; // Check

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
    //fprintf(stderr, "\nb dot b = %+.15e", inner_product(X_u, b_u, b_u));

    // 3. If spacelike, OK.
    // 4. If not, flip sign of b_u[0]? Or just create new b_u?
    while(inner_product(X_u, b_u, b_u) < 0.1)
    {
        LOOP_i b_u[i] = 1. / sqrt(3.) * genrand_RCARRY();
    }

    // Now we have b_u s.t. b dot b is spacelike. CHECK!


    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);

    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] = (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive; // CHECK

    // Permutation symbol eta, [ijkl]. Can be made into the contravariant Levi-Civita tensor by multiplying with -1/sqrt(g), g being the determinant of the covariant metric.
    // (source: cpluscode.blogspot.nl/p/physics-c.html)
    double g = determ(g_dd, 4);
    double eps[4][4][4][4];
    LOOP_ijkl{
        if((i==j) || (i==k) || (i==l) || (j==k) || (j==l) || (k==l))
            eps[i][j][k][l] = 0.;
        else eps[i][j][k][l] =
            ((i-j) * (i-k) * (i-l) * (j-k) * (j-l) * (k-l) / 12.);
    }

//    fprintf(stderr, "\n\nU DOT T: %+.15e\n\n", inner_product(X_u, U_u, b_u));
//    fprintf(stderr, "k DOT T: %+.15e\n\n", inner_product(X_u, k_u, b_u));
//    fprintf(stderr, "k DOT k: %+.15e\n\n", inner_product(X_u, k_u, k_u));
//    fprintf(stderr, "U DOT U: %+.15e\n\n", inner_product(X_u, U_u, U_u));

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i    e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] += (-1./sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) / (omega * Ncursive); // Check



    double r_current = logscale ? exp(X_u[1]) : X_u[1];
//    photon_u[5] = (logscale ? 1. / rPhoton * k_u1 : k_u1); // Covers BL and MBL


    // Construct the tetrad with contravariant coordinate index
    // CONVENTION: t, perp, para, K <=> t, x, y, z
    tetrad_u[0][0] = e_u_t[0];
    tetrad_u[1][0] = e_u_t[1];
    tetrad_u[2][0] = e_u_t[2];
    tetrad_u[3][0] = e_u_t[3];

    tetrad_u[0][2] = e_u_para[0];
    tetrad_u[1][2] = e_u_para[1];
    tetrad_u[2][2] = e_u_para[2];
    tetrad_u[3][2] = e_u_para[3];

    tetrad_u[0][1] = e_u_perp[0];
    tetrad_u[1][1] = e_u_perp[1];
    tetrad_u[2][1] = e_u_perp[2];
    tetrad_u[3][1] = e_u_perp[3];

    tetrad_u[0][3] = e_u_K[0];
    tetrad_u[1][3] = e_u_K[1];
    tetrad_u[2][3] = e_u_K[2];
    tetrad_u[3][3] = e_u_K[3];
}




// Here we implement Ziri's tetrad
void create_observer_tetrad_u2(const double X_u[], const double k_u[], const double U_u[], const double b_u[], double tetrad_u[][4])
{
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
    LOOP_i e_u_K[i] = k_u[i] / omega - U_u[i]; // Check

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

    // Only requirement: b dot b > 0 (i.e., b is spacelike) and b points along Y direction of image.

//    double b_u[4] = {0., 0., 1., 0.};

    // Now we have b_u s.t. b dot b is spacelike. CHECK!


    // First some required quantities:
    double Beta = inner_product(X_u, U_u, b_u);
    double Ccursive = inner_product(X_u, k_u, b_u) / omega - Beta;
    double b2 = inner_product(X_u, b_u, b_u);
    double Ncursive = sqrt(b2 + Beta * Beta - Ccursive * Ccursive);


    // Now we can construct e_u_para:
    LOOP_i e_u_para[i] = (b_u[i] + Beta * U_u[i] - Ccursive * e_u_K[i]) / Ncursive; // CHECK

    // Permutation symbol eta, [ijkl]. Can be made into the contravariant Levi-Civita tensor by multiplying with -1/sqrt(g), g being the determinant of the covariant metric.
    // (source: cpluscode.blogspot.nl/p/physics-c.html)
    double g = determ(g_dd, 4);
    double eps[4][4][4][4];
    LOOP_ijkl{
        if((i==j) || (i==k) || (i==l) || (j==k) || (j==l) || (k==l))
            eps[i][j][k][l] = 0.;
        else eps[i][j][k][l] =
            ((i-j) * (i-k) * (i-l) * (j-k) * (j-l) * (k-l) / 12.);
    }

//    fprintf(stderr, "\n\nU DOT T: %+.15e\n\n", inner_product(X_u, U_u, b_u));
//    fprintf(stderr, "k DOT T: %+.15e\n\n", inner_product(X_u, k_u, b_u));
//    fprintf(stderr, "k DOT k: %+.15e\n\n", inner_product(X_u, k_u, k_u));
//    fprintf(stderr, "U DOT U: %+.15e\n\n", inner_product(X_u, U_u, U_u));

    // Need b_d
    double b_d[4];
    lower_index(X_u, b_u, b_d);

    LOOP_i    e_u_perp[i] = 0.;
    LOOP_ijkl e_u_perp[i] += (-1./sqrt(-g) * eps[i][j][k][l] * U_d[j] * k_d[k] * b_d[l]) / (omega * Ncursive); // Check



    double r_current = 1;//logscale ? exp(X_u[1]) : X_u[1];


//    photon_u[5] = (logscale ? 1. / rPhoton * k_u1 : k_u1); // Covers BL and MBL

    // Construct the tetrad with contravariant coordinate index
    // CONVENTION: t, perp, para, K <=> t, x, y, z
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








double tetrad_identity_eta(const double X_u[4], const double tetrad_u[4][4], const int a, const int b){
    double result = 0.;

    double g_dd[4][4];
    metric_dd(X_u, g_dd);

    int i, j;
    LOOP_ij result += g_dd[i][j] * tetrad_u[i][a] * tetrad_u[j][b];

    return result;
}

double tetrad_identity_g(const double tetrad_u[][4], const int mu, const int nu){

    double eta_uu[4][4] = {{-1., 0., 0., 0.},
        {0., 1., 0., 0.},
        {0., 0., 1., 0.},
        {0., 0., 0., 1.}};

    double result = 0.;

    int i, j;
    LOOP_ij result += eta_uu[i][j] * tetrad_u[mu][i] * tetrad_u[nu][j];

    return result;
}

double tetrad_identity_sum_latin(const double tetrad_u[4][4], const double tetrad_d[4][4], const int mu, const int nu){

    double result = 0.;

    int i;
    LOOP_i result += tetrad_u[mu][i] * tetrad_d[nu][i];

    return result;
}

double tetrad_identity_sum_greek(const double tetrad_u[4][4], const double tetrad_d[4][4], const int a, const int b){

    double result = 0.;

    int i;
    LOOP_i result += tetrad_u[i][a] * tetrad_d[i][b];

    return result;
}

void create_tetrad_d(const double X_u[], const double tetrad_u[][4], double tetrad_d[][4]){
    double eta_minkowski[4][4] = {
        {-1., 0., 0., 0.},
        {0.,  1., 0., 0.},
        {0.,  0., 1., 0.},
        {0.,  0., 0., 1.},
    };

    int i,j,q,s, k, l;

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    // Create the tetrad with covariant coordinate index:
    LOOP_ij tetrad_d[i][j] = 0.;

    // ***NOTE*** the index order must be swapped here, i.e. we take the transpose.
    LOOP_ij{
        LOOP_kl tetrad_d[j][i] += eta_minkowski[i][k] * g_dd[j][l] * tetrad_u[l][k];
    }
}


void check_tetrad_identities(const double X_u[], double tetrad_u[][4]){
    printf("\nMinkowski");

    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 0) - (-1.)), fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 1)), fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 2)), fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 3)));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 0)), fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 1) - 1.), fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 2)), fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 3)));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 0)), fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 1)), fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 2) - 1.), fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 3)));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 0)), fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 1)), fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 2)), fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 3) - 1.));

    printf("\ng");

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double tetrad_d[4][4];
    create_tetrad_d(X_u, tetrad_u, tetrad_d);

    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_g(tetrad_u, 0, 0) - g_uu[0][0]), fabs(tetrad_identity_g(tetrad_u, 0, 1) - g_uu[0][1]), fabs(tetrad_identity_g( tetrad_u, 0, 2) - g_uu[0][2]), fabs(tetrad_identity_g( tetrad_u, 0, 3) - g_uu[0][3]));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_g(tetrad_u, 1, 0) - g_uu[1][0]), fabs(tetrad_identity_g(tetrad_u, 1, 1) - g_uu[1][1]), fabs(tetrad_identity_g( tetrad_u, 1, 2) - g_uu[1][2]), fabs(tetrad_identity_g( tetrad_u, 1, 3) - g_uu[1][3]));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_g(tetrad_u, 2, 0) - g_uu[2][0]), fabs(tetrad_identity_g(tetrad_u, 2, 1) - g_uu[2][1]), fabs(tetrad_identity_g( tetrad_u, 2, 2) - g_uu[2][2]), fabs(tetrad_identity_g( tetrad_u, 2, 3) - g_uu[2][3]));
    printf("\n%.3e %.3e %.3e %.3e", fabs(tetrad_identity_g( tetrad_u, 3, 0) - g_uu[3][0]), fabs(tetrad_identity_g( tetrad_u, 3, 1) - g_uu[3][1]), fabs(tetrad_identity_g( tetrad_u, 3, 2) - g_uu[3][2]), fabs(tetrad_identity_g( tetrad_u, 3, 3) - g_uu[3][3]));

    printf("\ndelta1 (sum latin)");

    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 0), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 1), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 2), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 0), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 1), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 2), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 0), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 1), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 2), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 0), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 1), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 2), tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 3));

    printf("\ndelta2 (sum greek)");

    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 0), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 1), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 2), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 0), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 1), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 2), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 0), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 1), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 2), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 3));
    printf("\n%+.15e %+.15e %+.15e %+.15e", tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 0), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 1), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 2), tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 3));

    printf("\n");
}


void check_tetrad_compact(const double X_u[], const double tetrad_u[][4])
{
    double result = 0.;

    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 0) - (-1.)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 1)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 2)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 0, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 0)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 1) - 1.) + fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 2)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 1, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 0)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 1)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 2) - 1.) + fabs(tetrad_identity_eta(X_u, tetrad_u, 2, 3));
    result += fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 0)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 1)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 2)) + fabs(tetrad_identity_eta(X_u, tetrad_u, 3, 3) - 1.);

    // Obtain relevant metric terms:
    double g_uu[4][4], g_dd[4][4];
    metric_uu(X_u, g_uu);
    metric_dd(X_u, g_dd);

    double tetrad_d[4][4];
    int i, j;
    LOOP_ij tetrad_d[i][j] = 0.;
    create_tetrad_d(X_u, tetrad_u, tetrad_d);

    result += fabs(tetrad_identity_g( tetrad_u, 0, 0) - g_uu[0][0]) + fabs(tetrad_identity_g( tetrad_u, 0, 1) - g_uu[0][1]) + fabs(tetrad_identity_g( tetrad_u, 0, 2) - g_uu[0][2]) + fabs(tetrad_identity_g( tetrad_u, 0, 3) - g_uu[0][3]);
    result += fabs(tetrad_identity_g( tetrad_u, 1, 0) - g_uu[1][0]) + fabs(tetrad_identity_g( tetrad_u, 1, 1) - g_uu[1][1]) + fabs(tetrad_identity_g( tetrad_u, 1, 2) - g_uu[1][2]) + fabs(tetrad_identity_g( tetrad_u, 1, 3) - g_uu[1][3]);
    result += fabs(tetrad_identity_g( tetrad_u, 2, 0) - g_uu[2][0]) + fabs(tetrad_identity_g( tetrad_u, 2, 1) - g_uu[2][1]) + fabs(tetrad_identity_g( tetrad_u, 2, 2) - g_uu[2][2]) + fabs(tetrad_identity_g( tetrad_u, 2, 3) - g_uu[2][3]);
    result += fabs(tetrad_identity_g( tetrad_u, 3, 0) - g_uu[3][0]) + fabs(tetrad_identity_g( tetrad_u, 3, 1) - g_uu[3][1]) + fabs(tetrad_identity_g( tetrad_u, 3, 2) - g_uu[3][2]) + fabs(tetrad_identity_g( tetrad_u, 3, 3) - g_uu[3][3]);

    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 0) - 1.) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 1)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 2)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 0, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 0)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 1) - 1.) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 2)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 1, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 0)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 1)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 2) - 1.) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 2, 3));
    result += fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 0)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 1)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 2)) + fabs(tetrad_identity_sum_greek(tetrad_u, tetrad_d, 3, 3) - 1.);

    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 0) - 1.) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 1)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 2)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 0, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 0)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 1) - 1.) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 2)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 1, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 0)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 1)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 2) - 1.) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 2, 3));
    result += fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 0)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 1)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 2)) + fabs(tetrad_identity_sum_latin(tetrad_u, tetrad_d, 3, 3) - 1.);

    printf("\nTetrad identities (should be close to zero): %+.15e", result);
}
















// Return flux as a function of radius
double no_torque_flux(double X_u[4])
{
    double r_current = logscale ? exp(X_u[1]) : X_u[1];

    double G_gravity = 1.;
    double M_mass = 1.;
    double Mdot_0 = 1.; // Accretion rate

    double f_function = 1.;

    double r3 = r_current * r_current * r_current;

    double x0 = sqrt(R_ISCO);
    double x = sqrt(r_current);
    double x1 = 2. * cos(1. / 3. * acos(a) - M_PI / 3.); // CHANCE OF MISTAKE - the cos^-1
    double x2 = 2. * cos(1. / 3. * acos(a) + M_PI / 3.); // CHANCE OF MISTAKE - the cos^-1
    double x3 = -2. * cos(1. / 3. * acos(a)); // CHANCE OF MISTAKE - the cos^-1

    double B_cursive = 1. + a * pow(x, -3.);
    double C_cursive = 1. - 3. * pow(x, -2.) + 2. * a * pow(x, -3.);
    double Q_cursive = (1. + a * pow(x, -3.)) / (x * (1. - 3. * pow(x, -2) + 2. * a * pow(x, -3.))) * // Initial fraction
                       (x - x0 - 1.5 * a * log(x / x0) - // Square-bracket term starts on this line
                       3. * (x1 - a) * (x1 - a) / (x1 * (x1 - x2) * (x1 - x3)) * log((x - x1) / (x0 - x1)) -
                       3. * (x2 - a) * (x2 - a) / (x2 * (x2 - x1) * (x2 - x3)) * log((x - x2) / (x0 - x2)) -
                       3. * (x3 - a) * (x3 - a) / (x3 * (x3 - x1) * (x3 - x2)) * log((x - x3) / (x0 - x3)));

    double F = 3. * G_gravity * M_mass * Mdot_0 * f_function / (8. * M_PI * r3) * Q_cursive * B_cursive / sqrt(C_cursive);

    return F;
}

double no_torque_flux_ben(double X_u[4]){
    double r_current = logscale ? exp(X_u[1]) : X_u[1];
    r_current;
    double x = sqrt(r_current);
    double x0 = sqrt(R_ISCO);
    double x1 = 2. * cos(1. / 3. * acos(a) - M_PI / 3.);
    double x2 = 2. * cos(1. / 3. * acos(a) + M_PI / 3.);
    double x3 = -2. * cos(1. / 3. * acos(a));

    double f = x * x / (r_current * r_current * r_current) * 1. / (x * x * x - 3. * x + 2. * a) * 
               (x - x0 - 1.5 * a * log(x / x0) - // Square-bracket term starts on this line
                       3. * (x1 - a) * (x1 - a) / (x1 * (x1 - x2) * (x1 - x3)) * log((x - x1) / (x0 - x1)) -
                       3. * (x2 - a) * (x2 - a) / (x2 * (x2 - x1) * (x2 - x3)) * log((x - x2) / (x0 - x2)) -
                       3. * (x3 - a) * (x3 - a) / (x3 * (x3 - x1) * (x3 - x2)) * log((x - x3) / (x0 - x3)));

//    return 6.87981e24 * f;
    return 3. / (8. * M_PI) * GGRAV * MBH * 1.399e17 * pow(GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT), -3.) * f;
}





// Get params from Chandrasekhar (1960)
void get_IF_and_p(const double mu, double *I_over_F, double *deg_of_p)
{

/*
mu        I/F            p
0.00    0.41441        0.11713
0.05    0.47490        0.08979
0.10    0.52397        0.07448
0.15    0.57001        0.06311
0.20    0.61439        0.05410
0.25    0.65770        0.04667
0.30    0.70029        0.04041
0.35    0.74234        0.03502
0.40    0.78398        0.03033
0.45    0.82530        0.02619
0.50    0.86637        0.02252
0.55    0.90722        0.01923
0.60    0.94789        0.01627
0.65    0.98842        0.01358
0.70    1.02882        0.011123
0.75    1.06911        0.008880
0.80    1.10931        0.006818
0.85    1.14943        0.004919
0.90    1.18947        0.003155
0.95    1.22945        0.001522
1.00    1.26938        0.0
*/

    double IF_values[21] = {
0.41441,
0.47490,
0.52397,
0.57001,
0.61439,
0.65770,
0.70029,
0.74234,
0.78398,
0.82530,
0.86637,
0.90722,
0.94789,
0.98842,
1.02882,
1.06911,
1.10931,
1.14943,
1.18947,
1.22945,
1.26938};

    double p_values[21] = {
0.11713,
0.08979,
0.07448,
0.06311,
0.05410,
0.04667,
0.04041,
0.03502,
0.03033,
0.02619,
0.02252,
0.01923,
0.01627,
0.01358,
0.011123,
0.008880,
0.006818,
0.004919,
0.003155,
0.001522,
0.0};


    // Nearest neighbor interpolation
    // Convert from mu to index
    int index = (int) (mu * 20.);

    *I_over_F = IF_values[index];
    *deg_of_p = p_values[index];

    // Linear interpolation
    double mu_min = 0.;
    double mu_max = 1.;
    double Dmu = 0.05;

    int lower_index = floor(mu / Dmu);
    int higher_index = lower_index + 1; // TODO: TO DO: WHAT IF  lower index + 1 exceeds array....

    double dx_interp = (mu - (double) lower_index * Dmu) / Dmu;

    *I_over_F = IF_values[lower_index] * (1. - dx_interp) + IF_values[higher_index] * dx_interp;
    *deg_of_p = p_values[lower_index] * (1. - dx_interp) + p_values[higher_index] * dx_interp;
}







// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_parallel(const double y[], double complex f_u[], double fvector[], double complex f_u_vector[]){
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];
    int i, j, k; // Einstein summation over indices v and w

    LOOP_ijk gamma_udd[i][j][k] = 0.;

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double complex A_u[4] = {0.,   0.,   0.,   0.  }; // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
    connection_udd(X_u,gamma_udd);
    //connection_num_udd(X_u, gamma_udd);

    // Compute 4-acceleration using the geodesic equation
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];
    LOOP_i {
        fvector[i]     = U_u[i];
        fvector[i + 4] = A_u[i];
    }

    // Reset A_u
    LOOP_i A_u[i] = 0.;

    // Compute f_u vector acceleration
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * f_u[k];
    LOOP_i {
        f_u_vector[i] = A_u[i];
    }
}

void rk4_step_f(double y[], double complex f_u[], double dt){
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
    double complex weights[4] = {0.5, 0.5, 1., 0.}; // Weights used for updating y
    for (q = 0; q < 4; q++) {
        f_parallel(yshift, f_u_shift, fvector, f_u_vector); // Apply function f to current y to obtain fvector
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
    for (i = 0; i < 4 * 2; i++){
        y[i] = y[i] + 1. / 6. * (dx[0 * 4 * 2 + i] +
                                 dx[1 * 4 * 2 + i] * 2. +
                                 dx[2 * 4 * 2 + i] * 2. +
                                 dx[3 * 4 * 2 + i]);
    }

    // Update the f-vector (polarization)
    for (i = 0; i < 4; i++){
        f_u[i] = f_u[i] + 1. / 6. * (df[0 * 4 + i] +
                                     df[1 * 4 + i] * 2. +
                                     df[2 * 4 + i] * 2. +
                                     df[3 * 4 + i]);
    }
}


/*
double complex inner_product_complex_complex(const double *X_u, double complex *A_u, double complex *B_u){
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

double complex inner_product_real_complex(const double *X_u, double *A_u, double complex *B_u){
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
*/
// NOTE: works only in Kerr metric
// Ziri's suggestion: construct U vecs
void construct_U_vector(const double X_u[], double U_u[])
{
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
    //printf("\nU dot U: %+.15e", inner_product(X_u, U_u, U_u));
}





void f_tetrad_u_to_stokes(const double complex f_tetrad_u[4], const double p, double complex S_A[4]){
    double complex IfromJones = 0.;
    double complex QfromJones = 0.;//cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) - cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]);
    double complex UfromJones = 0.;//conj(f_tetrad_u[1]) * f_tetrad_u[2] + f_tetrad_u[1] * conj(f_tetrad_u[2]);
    double complex VfromJones = 0.;//I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] - f_tetrad_u[1] * conj(f_tetrad_u[2]));

    IfromJones = (cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) + cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2])) / p;
    QfromJones = cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) - cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]);
    UfromJones = conj(f_tetrad_u[1]) * f_tetrad_u[2] + f_tetrad_u[1] * conj(f_tetrad_u[2]);
    VfromJones = I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] - f_tetrad_u[1] * conj(f_tetrad_u[2]));



    S_A[0] = IfromJones;
    S_A[1] = QfromJones;
    S_A[2] = UfromJones;
    S_A[3] = VfromJones;

    if(S_A[0] != S_A[0] || S_A[1] != S_A[1] ||S_A[2] != S_A[2] ||S_A[3] != S_A[3]){
        printf("\n AND P IS %+.15e", p);
        printf("\nf_u[1] = %+.15e, %+.15e; f_u[2] = %+.15e, %+.15e", creal(f_tetrad_u[1]), cimag(f_tetrad_u[1]), creal(f_tetrad_u[2]), cimag(f_tetrad_u[2]));
        printf("NANS IN DA HOUSE \n NANS IN DA HOUSE \n NANS IN DA HOUSE");
    }
}








// Integrate the null geodesic defined by "photon_u"
double integrate_geodesic(double alpha, double beta, double *photon_u, double *lightpath, int *steps, double cutoff_inner, double *f_x, double *f_y, double *p, double *IQUV){
    int i, q;
    double t_init = 0.;
    double dlambda_adaptive;
    int theta_turns = 0;
    double thetadot_prev = 0.;
    double X_u[4], k_u[4];

    // Create initial ray conditions
    initialize_photon_parallel(alpha, beta, photon_u, t_init);
  //  initialize_photon(alpha, beta, photon_u, t_init);

    // Current r-coordinate
    double r_current = logscale ? exp(photon_u[1]) : photon_u[1];

    // Reset lambda and steps
    double lambda = 0.;
    *steps = 0;

    double theta_prev = 0.;
    int TERMINATE = 0; // Termination condition for ray

    // Trace light ray until it reaches the event horizon or the outer
    // cutoff, or steps > max_steps
    // Stop condition for BL coords
    while (r_current > cutoff_inner && r_current < cutoff_outer &&
           *steps < max_steps && !TERMINATE)
    {
        // Current photon position/wave vector
        LOOP_i{
            X_u[i] = photon_u[i];
            k_u[i] = photon_u[i + 4];
        }

        // Possibly terminate ray to eliminate higher order images
        if (thetadot_prev * photon_u[6] < 0. && *steps > 2){
            theta_turns += 1;
        }
        thetadot_prev = photon_u[6];
        if((beta < 0. && theta_turns > max_order) || (beta > 0. && theta_turns > (max_order + 1)))
            TERMINATE = 1;

        // Compute adaptive step size
        //dlambda_adaptive = -STEPSIZE;
        dlambda_adaptive = stepsize(X_u, k_u);

        // Enter current position/velocity/dlambda into lightpath
        for (q = 0; q < 8; q++)
            lightpath[*steps * 9 + q] = photon_u[q];
        lightpath[*steps * 9 + 8] = fabs(dlambda_adaptive);

        // Advance ray/particle
	#if(int_method == RK4)
            rk4_step(photon_u, &f_geodesic, dlambda_adaptive);
	#elif(int_method == VER)
            verlet_step(photon_u, &f_geodesic, dlambda_adaptive);
	#endif

        // Advance (affine) parameter lambda and 'steps' variable
        lambda += fabs(dlambda_adaptive);
        r_current = logscale ? exp(photon_u[1]) : photon_u[1] ;

        *f_x = 0.;
        *f_y = 0.;

        // THIN-DISK MODEL
        //////////////////

	// Get the theta coordinate right for MKS2 (where theta is between 0 and 1) and everything else
	double factor = 1.;
        #if(metric == MKS2)
	    factor = M_PI;
        #endif

	// Disk inner and outer radius. For the N-K model, the outer radius is large.
	double R_inner = R_ISCO;
        double R_outer = 100.;

        int i,j;

        double tetrad_u[4][4], tetrad_d[4][4];
        LOOP_ij tetrad_u[i][j] = 0.;
        LOOP_ij tetrad_d[i][j] = 0.;

        double k_u_tetrad[4], n_u[4], n_u_tetrad[4];
        LOOP_i {
            n_u[i] = 0.;
            k_u_tetrad[i] = 0.;
            n_u_tetrad[i] = 0.;
        }
        n_u[2] = -1.; // Normal vector for disk. Points in theta direction

        double k_cart[3] = {0., 0., 0.};
        double n_cart[3] = {0., 0., 0.};
        double complex f_cart[3] = {0., 0., 0.};
        double complex f_u_tetrad[4] = {0., 0., 0., 0.};
        double complex f_u[4] = {0., 0., 0., 0.};

        if((photon_u[2] * factor - M_PI / 2.) * (theta_prev * factor - M_PI / 2.) < 0. && (r_current < R_inner || r_current > R_outer))
            return 0.;

        // If we are in the disk...
        if((photon_u[2] * factor - M_PI / 2.) * (theta_prev * factor - M_PI / 2.) < 0. &&
           r_current > R_inner && r_current < R_outer)
        {
            double j_nu, Uplasma_u[4], k_d[4];
            double nu_p;
            LOOP_i k_d[i] = 0.;
            LOOP_i Uplasma_u[i] = 0.;

            LOOP_i {
                k_u_tetrad[i] = 0.;
                n_u_tetrad[i] = 0.;
            }
            LOOP_ij tetrad_u[i][j] = 0.;
            LOOP_ij tetrad_d[i][j] = 0.;


            double THINDISK_TEST_FREQ = 2.417989e17;

            // Scale the wave vector
            LOOP_i k_u[i] *= PLANCK_CONSTANT * THINDISK_TEST_FREQ /
                             (ELECTRON_MASS * SPEED_OF_LIGHT*SPEED_OF_LIGHT);

            //lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Obtain disk velocity
            disk_velocity(X_u, Uplasma_u);

            // SCALE VELOCITY VECTOR?!?!?!?!?!!??!?!?!
            double Uplasma_u_scaled[4] = {0., 0., 0., 0.};
            LOOP_i Uplasma_u_scaled[i] = Uplasma_u[i];

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(Uplasma_u_scaled, k_d);//_natural_units(Uplasma_u, k_d);


	    // THIN DISK EMISSION
            /////////////////////

            double Fcurs = no_torque_flux_ben(X_u);
            double Teff = pow(Fcurs / STEFAN_BOLTZMANN, 1./4.);
            double n_hard = 1.8;
            double I_nu_local = 1. / (pow(n_hard, 4)) * planck_function(nu_p, n_hard * Teff);

            // POLARIZATION VECTOR
            //////////////////////

//            printf("k dot k = %+.15e\n", inner_product(X_u, k_u, k_u));

            create_tetrad_u2(X_u, k_u, Uplasma_u, tetrad_u);
            //check_tetrad_compact(X_u, tetrad_u);

            create_tetrad_d(X_u, tetrad_u, tetrad_d);

            // Compute k_u_tetrad and n_u_tetrad
            LOOP_ij{
                k_u_tetrad[i] += tetrad_d[j][i] * k_u[j];
                n_u_tetrad[i] += tetrad_d[j][i] * n_u[j];
            }

//	    printf("\nk_u_tetrad norm = %+.15e", -k_u_tetrad[0] * k_u_tetrad[0] + k_u_tetrad[1] * k_u_tetrad[1] + k_u_tetrad[2] * k_u_tetrad[2] + k_u_tetrad[3] * k_u_tetrad[3]);

            k_cart[0] = k_u_tetrad[1];
            k_cart[1] = k_u_tetrad[2];
            k_cart[2] = k_u_tetrad[3];

            n_cart[0] = n_u_tetrad[1];
            n_cart[1] = n_u_tetrad[2];
            n_cart[2] = n_u_tetrad[3];

            double k_norm = sqrt(k_cart[0] * k_cart[0] + k_cart[1] * k_cart[1] + k_cart[2] * k_cart[2]);
            double n_norm = sqrt(n_cart[0] * n_cart[0] + n_cart[1] * n_cart[1] + n_cart[2] * n_cart[2]);

            LOOP_i if(i < 3) k_cart[i] /= k_norm;
            LOOP_i if(i < 3) n_cart[i] /= n_norm;

            double k_dot_n = k_cart[0] * n_cart[0] + k_cart[1] * n_cart[1] + k_cart[2] * n_cart[2];

            // Normal vector should always point along with wave vector. Flip if it is not the case.
            if(k_dot_n < 0.){
                k_dot_n *= -1.;
                n_cart[0] *= -1.;
                n_cart[1] *= -1.;
                n_cart[2] *= -1.;
                n_u[0] *= -1.;
                n_u[1] *= -1.;
                n_u[2] *= -1.;
                n_u[3] *= -1.;
            }

            f_cart[0] = (k_cart[1] * n_cart[2] - k_cart[2] * n_cart[1]);
            f_cart[1] = (k_cart[2] * n_cart[0] - k_cart[0] * n_cart[2]);
            f_cart[2] = (k_cart[0] * n_cart[1] - k_cart[1] * n_cart[0]);

            // NORMALIZE the polarization vector (cross product of two normal vectors is not always normalized)
            double f_norm = sqrt(f_cart[0] * f_cart[0] + f_cart[1] * f_cart[1] + f_cart[2] * f_cart[2]);
            LOOP_i if(i < 3) f_cart[i] /= f_norm;

            f_u_tetrad[0] = 0.;
            f_u_tetrad[1] = f_cart[0];
            f_u_tetrad[2] = f_cart[1];
            f_u_tetrad[3] = 0.;

            // Compute f_u
            LOOP_i f_u[i] = 0.;
            LOOP_ij f_u[i] += tetrad_u[i][j] * f_u_tetrad[j];

            // Now we have f_u_i, the initial polarization vector at the disk.

            // Must obtain p (degree of polarization) and I/F (limb-darkening) from Chandrasekhar tables...
            double mu = k_dot_n;
            double deg_of_p = 0.;
            double I_over_F = 0.;
            get_IF_and_p(k_dot_n, &I_over_F, &deg_of_p);

            double I_inv = I_over_F * I_nu_local / nu_p / nu_p / nu_p;
            double I_at_cam = I_inv * THINDISK_TEST_FREQ * THINDISK_TEST_FREQ * THINDISK_TEST_FREQ;


            // POLARIZED TRANSFER
            /////////////////////

            // If p > minimum...

            // Parallel-transport f_u along lightpath to X_u_0.
            int index;
            double X_u_current[4] = {0., 0., 0., 0.};
            double k_u_current[4] = {0., 0., 0., 0.};
            double photon_u_current[8] = {0., 0., 0., 0., 0., 0., 0., 0.};

            double dl_current = 0.;

            // Perform backward parallel transport of the polarization vector.
            for(index = *steps; index > -1; index--)
            {
                for (q = 0; q < 8; q++){
                    photon_u_current[q] = lightpath[index * 9 + q];
                }
                dl_current = lightpath[(index - 1) * 9 + 8];

                if(index > 0)
                    rk4_step_f(photon_u_current, f_u, dl_current);
            }

            LOOP_i X_u_current[i] = photon_u_current[i];
            LOOP_i k_u_current[i] = photon_u_current[i + 4];

            // Construct the observer tetrad.
            // X_u_current and k_u_current are simply the initial position and wave vector.
            // Note that k_u_current points INTO the camera sensor plane.
            double cam_up_u[4] = {0., 0., 0., -1.};

            // Need U_obs_u
            double U_obs_u[4] = {0., 0., 0., 0.};
            double obs_tetrad_u[4][4], obs_tetrad_d[4][4];
            LOOP_ij obs_tetrad_u[i][j] = 0.;
            LOOP_ij obs_tetrad_d[i][j] = 0.;
            construct_U_vector(X_u_current, U_obs_u);


            double k_u_current_reverse[4] = {k_u_current[0], k_u_current[1], k_u_current[2], k_u_current[3]};
            k_u_current_reverse[1] *= -1.;
            k_u_current_reverse[2] *= -1.;
            k_u_current_reverse[3] *= -1.;

            normalize_null(X_u_current, k_u_current_reverse);

            create_observer_tetrad_u2(X_u_current, k_u_current, U_obs_u, cam_up_u, obs_tetrad_u);


            create_tetrad_d(X_u_current, obs_tetrad_u, obs_tetrad_d);

            // Compute the final polarization state.
            double complex f_u_obs_tetrad[4] = {0., 0., 0., 0.};
            LOOP_ij f_u_obs_tetrad[i] += obs_tetrad_d[j][i] * creal(f_u[j]);

            *p = deg_of_p;
            *f_x = f_u_obs_tetrad[1];
            *f_y = f_u_obs_tetrad[2];

            double complex S_A[4] = {0., 0., 0., 0.};


            f_tetrad_u_to_stokes(f_u_obs_tetrad, 1., S_A);

            double I_tot = I_at_cam;//I_inv * I_over_F;
            double I_pol = I_tot * *p;

            // Construct final (NON-INVARIANT) Stokes params.
            double S_If = I_tot;
            double S_Qf = S_A[1] * I_pol;
            double S_Uf = S_A[2] * I_pol;
            double S_Vf = S_A[3] * I_pol;

            IQUV[0] = S_If;
            IQUV[1] = S_Qf;
            IQUV[2] = S_Uf;
            IQUV[3] = S_Vf;







            // Return intensity.
            return 0.;//I_inv * I_over_F;
        }

//        if((photon_u[2] * factor - M_PI / 2.) * (theta_prev * factor - M_PI / 2.) < 0.)
  //          return 0.;

        // Update steps variable.
        *steps = *steps + 1;

 	// Update 'theta_prev' to check whether we cross the equatorial plane in this step.
	theta_prev = photon_u[2];
    }

    IQUV[0] = 0.;
    IQUV[1] = 0.;
    IQUV[2] = 0.;
    IQUV[3] = 0.;

    return 0.;
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

    double a_nu = 0.;

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
          // pitch_ang = 3.1415/3.;

            // CGS UNITS USED FROM HERE ON OUT
            //////////////////////////////////

            // Scale the wave vector to correct energy
            LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                             (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // Convert distance dlambda accordingly
            dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) / (PLANCK_CONSTANT * frequency);

            //lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            //nu_p = CAM_FREQ;
            nu_p = freq_in_plasma_frame(Uplasma_u, k_d);
            nu_p2 = nu_p * nu_p;

            // Obtain emission coefficient in current plasma conditions
           j_nu = emission_coeff_THSYNCHAV(B, THETA_e, nu_p, n_e);

            // Obtain absorption coefficient
            if (ABSORPTION){
                a_nu = absorption_coeff_TH(j_nu, nu_p, THETA_e);
            }

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT / (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // Solve the transport equation with emission and self-absoprtion
            dI = (j_nu / nu_p2 - nu_p * a_nu * I_current) * dl_current * C;


            // Only add I_current if it is not NaN
            if(j_nu == j_nu && exp(X_u[1]) < RT_OUTER_CUTOFF){ // I_current += exp(-tau) * j_nu / nu_p / nu_p * dl_current * C;
                I_current += dI; // Old way of integrating
             }
        }
    }

    // Store integrated intensity in the image
    return I_current * pow(frequency, 3.);
}
*/
