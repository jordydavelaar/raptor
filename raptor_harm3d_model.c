/*
 * raptor_harm_model.c
 *
 * Please note that most of the code in this file was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 */

#include "raptor_harm3d_model.h"
#include "constants.h"
#include "functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameters.h"

/* HDF5 v1.8 API */

void init_model() {
    /* find dimensional quantities from black hole
     mass and its accretion rate */
    set_units(M_UNIT);

    fprintf(stderr, "getting simulation data...\n");

    init_harm3d_data(GRMHD_FILE);
}

void init_harm3d_data(char *fname) {
    FILE *fp;

    fprintf(stderr, "%s\n", fname);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "can't open sim data file\n");
        exit(1);
    } else {
        fprintf(stderr, "successfully opened %s\n", fname);
    }

    /* get standard HARM header */
    fscanf(fp, "%d ", &N1);
    fscanf(fp, "%d ", &N2);
    fscanf(fp, "%d ", &N3);
    fscanf(fp, "%lf ", &gam);
    fscanf(fp, "%lf ", &a);
    fscanf(fp, "%lf ", &dx[1]);
    fscanf(fp, "%lf ", &dx[2]);
    fscanf(fp, "%lf ", &dx[3]);
    fscanf(fp, "%lf ", &startx[1]);
    fscanf(fp, "%lf ", &startx[2]);
    fscanf(fp, "%lf ", &startx[3]);
    fscanf(fp, "%lf ", &hslope);
    fscanf(fp, "%lf ", &Rin);
    fscanf(fp, "%lf ", &Rout);

    stopx[0] = 1.;
    stopx[1] = startx[1] + N1 * dx[1];
    stopx[2] = startx[2] + N2 * dx[2];
    stopx[3] = startx[3] + N3 * dx[3];
    fprintf(stderr, "phi limits is %e %e\n", startx[3], stopx[3]);

    init_storage();

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            for (int k = 0; k < N3; k++) {
                fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &p[KRHO][i][j][k],
                       &p[UU][i][j][k], &p[U1][i][j][k], &p[U2][i][j][k],
                       &p[U3][i][j][k], &p[B1][i][j][k], &p[B2][i][j][k],
                       &p[B3][i][j][k]);
            }
        }
    }

    fprintf(stderr, "%e\n", p[KRHO][N1 - 1][N2 - 1][N3 - 1]);

    fprintf(stderr, "done reading!\n");

    //	exit(1);
}


// Current metric: modified Kerr-Schild, squashed in theta
// to give higher resolution at the equator

// Get the fluid parameters - IN THE PLASMA FRAME?
void get_fluid_params(double X[NDIM], double *Ne, double *Thetae, double *B,
                      double Bcon[NDIM], double Ucon[NDIM], int *IN_VOLUME) {

    int i, j, k;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
    double gcon[NDIM][NDIM], gcov[NDIM][NDIM], Bcov[NDIM], Ucov[NDIM], coeff[4];
    double bsq, beta, beta_trans, b2, trat, two_temp_gam, Th_unit, Be;
    //double Rlow = 1, Rhigh = 1;
    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
        X[2] > stopx[2]) {
        *Ne = 0.;
        *IN_VOLUME = 0;
        return;
    }
    *IN_VOLUME = 1;

    Xtoijk(X, &i, &j, &k, del);

    metric_uu(X, gcon);
    metric_dd(X, gcov);

    coeff[1] = del[1];
    coeff[2] = del[2];
    coeff[3] = del[3];

    // now interpolated to geodesic location
    // interpolate (cubiclinear interp.) like in mibothros
    rho = interp_scalar(p[KRHO], i, j, k, coeff);
    uu = interp_scalar(p[UU], i, j, k, coeff);
    *Ne = rho * Ne_unit + 1e-40;

    // here unlike in mibothros it was interpolating scalars and
    // reconstructing velocity and magnetic field based on interpolated
    // coefficients
    Bp[1] = interp_scalar(p[B1], i, j, k, coeff);
    Bp[2] = interp_scalar(p[B2], i, j, k, coeff);
    Bp[3] = interp_scalar(p[B3], i, j, k, coeff);

    Vcon[1] = interp_scalar(p[U1], i, j, k, coeff);
    Vcon[2] = interp_scalar(p[U2], i, j, k, coeff);
    Vcon[3] = interp_scalar(p[U3], i, j, k, coeff);

    // reconstrueren van de 4 vectoren
    // Get Ucov
    VdotV = 0.;
    for (i = 1; i < NDIM; i++)
        for (j = 1; j < NDIM; j++)
            VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
    Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
    Ucon[0] = -Vfac * gcon[0][0];
    for (i = 1; i < NDIM; i++)
        Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

    lower_index(X, Ucon, Ucov);

    double Utot = 0;
    for (int i = 0; i < NDIM; i++)
        Utot += Ucon[i] * Ucov[i];

    // Get B and Bcov
    UdotBp = 0.;
    for (i = 1; i < NDIM; i++)
        UdotBp += Ucov[i] * Bp[i];
    Bcon[0] = UdotBp;
    for (i = 1; i < NDIM; i++)
        Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];

    lower_index(X, Bcon, Bcov);

    bsq = Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
          Bcon[3] * Bcov[3];

    *B = sqrt(bsq) * B_unit + 1e-40;

    /*electron temperature depending on the plasma magnetization*/
    beta = uu * (gam - 1.) / 0.5 / bsq;
    Be = (-(1. + gam * uu / (rho)) * Ucov[0]);
    beta_trans = 1.;
    b2 = pow(beta / beta_trans, 2);
    trat = 3.; // Rhigh * b2/(1. + b2) + Rlow /(1. + b2);
    two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
    Th_unit = (1.4444444444 - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + trat);
    *Thetae = (2. / 15.) * (uu / rho) * (PROTON_MASS / ELECTRON_MASS) + 1e-40;
    Be = (-(1. + two_temp_gam * uu / rho) * Ucov[0]);
    // if(bsq/rho>0.15){
    if (uu < 0)
        fprintf(stderr, "U %e %e\n", uu, p[UU][i][j][k]);
    ;

    if (*Thetae < 0)
        fprintf(stderr, "Te %e\n", *Thetae);
    if (*B < 0)
        fprintf(stderr, "B %e\n", *B);
    if (*Ne < 0)
        fprintf(stderr, "Ne %e %e\n", *Ne, p[KRHO][i][j][k]);

    if (bsq / rho > 1. || exp(X[1]) > 50.) {
        *Ne = 0;
        *IN_VOLUME = 0;
    }
}

void set_units(double M_unit_) {
    //	double MBH;

    /* set black hole mass */
    /** could be read in from file here,
        along with M_unit and other parameters **/
    //	MBH = 4.e6;

    /** input parameters appropriate to Sgr A* **/
    // double BH_MASS = MBH * MSUN;

    /** from this, calculate units of length, time, mass,
        and derivative units **/
    L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    T_unit = L_unit / SPEED_OF_LIGHT;

    fprintf(stderr, "\nUNITS\n");
    fprintf(stderr, "L,T,M: %g %g %g\n", L_unit, T_unit, M_unit_);

    RHO_unit = M_unit_ / pow(L_unit, 3);
    U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);

    fprintf(stderr, "rho,u,B: %g %g %g\n", RHO_unit, U_unit, B_unit);

    Ne_unit = RHO_unit / (PROTON_MASS + ELECTRON_MASS);
}

double interp_scalar(double ***var, int i, int j, int k, double coeff[4]) {

    double interp;
    int ip1, jp1, kp1;
    double b1, b2, b3, del[NDIM];

    del[1] = coeff[1];
    del[2] = coeff[2];
    del[3] = coeff[3];
    if (del[1] > 1 || del[2] > 1 || del[3] > 1 || del[1] < 0 || del[2] < 0 ||
        del[3] < 0)
        fprintf(stderr, "del[1] %e \n del[2] %e\n del[3] %e\n", del[1], del[2],
                del[3]);

    ip1 = i + 1;
    jp1 = j + 1;
    kp1 = k + 1;

    b1 = 1. - del[1];
    b2 = 1. - del[2];
    b3 = 1. - del[3];

    interp = var[i][j][k] * b1 * b2 + var[i][jp1][k] * b1 * del[2] +
             var[ip1][j][k] * del[1] * b2 + var[ip1][jp1][k] * del[1] * del[2];

    /* Now interpolate above in x3 */
    interp = b3 * interp + del[3] * (var[i][j][kp1] * b1 * b2 +
                                     var[i][jp1][kp1] * b1 * del[2] +
                                     var[ip1][j][kp1] * del[1] * b2 +
                                     var[ip1][jp1][kp1] * del[1] * del[2]);

    return interp;
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) {
    double phi;
    /* Map X[3] into sim range, assume startx[3] = 0 */
    phi = fmod(X[3], stopx[3]);
    if (phi < 0.)
        phi = stopx[3] + phi;

    *i = (int)((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
    *j = (int)((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
    *k = (int)(phi / dx[3] + 1000 - 0.5) - 1000;

    if (*i < 0) {
        *i = 0;
        del[1] = 0.;
    } else if (*i > N1 - 2) {
        *i = N1 - 2;
        del[1] = 1.;
    } else {
        del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
    }

    if (*j < 0) {
        *j = 0;
        del[2] = 0.;
    } else if (*j > N2 - 2) {
        *j = N2 - 2;
        del[2] = 1.;
    } else {
        del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
    }

    if (*k < 0) {
        *k = 0;
        del[3] = 0.;
    } else if (*k > N3 - 2) {
        *k = N3 - 2;
        del[3] = 1.;
    } else {
        del[3] = (phi - ((*k + 0.5) * dx[3])) / dx[3];
    }
    if (del[3] < 0)
        fprintf(stderr, "%e %e %e %d\n", del[3], phi, (*k + 0.5) * dx[3], *k);

    return;
}

// ALLOCATION STUFF BELOW HERE
//////////////////////////////

static void *malloc_rank1(int n1, int alloc_size) {
    void *A;

    if ((A = (void *)malloc(n1 * alloc_size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3) {
    double ***A;
    double *space;
    int i, j;

    space = malloc_rank1(n1 * n2 * n3, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = &(space[n3 * (j + n2 * i)]);
        }
    }

    return A;
}

void init_storage(void) {
    int i, j, k;

    p = (double ****)malloc(NPRIM * sizeof(double ***));
    for (i = 0; i < NPRIM; i++) {
        p[i] = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            p[i][j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                p[i][j][k] = (double *)malloc(N3 * sizeof(double));
            }
        }
    }

    fprintf(stderr, "done here with memory, %d %d %d %d\n", N1, N2, N3, NPRIM);
    return;
}
