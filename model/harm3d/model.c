/*
 * model file for HARM3D data
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////
double ****p;

int N1, N2, N3;

double gam;

double R0, Rin, Rout, a, hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];

double L_unit, T_unit;
double RHO_unit, U_unit, B_unit;
double Ne_unit, Thetae_unit;

// FUNCTIONS
////////////

void init_model() {

    set_units(M_UNIT);

    fprintf(stderr, "\nStarting read in of HARM3D GRMHD data...\n");

    init_grmhd_data(GRMHD_FILE);
}

void init_grmhd_data(char *fname) {
    FILE *fp;

    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
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
        if (i % (N1 / 3) == 0)
            fprintf(stderr, ".");
    }

    fprintf(stderr, "Done!\n");
}

int get_fluid_params(double X[NDIM], struct GRMHD *modvar) {

    int i, j, k;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], V_u[NDIM], Vfac, VdotV, UdotBp;
    double g_uu[NDIM][NDIM], g_dd[NDIM][NDIM], coeff[4];
    double bsq, beta, beta_trans, b2, trat, Th_unit, two_temp_gam;

    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
        X[2] > stopx[2]) {
        (*modvar).n_e = 0.;
        return 0;
    }

    Xtoijk(X, &i, &j, &k, del);

    metric_uu(X, g_uu);
    metric_dd(X, g_dd);

    coeff[1] = del[1];
    coeff[2] = del[2];
    coeff[3] = del[3];

    rho = interp_scalar(p[KRHO], i, j, k, coeff);
    uu = interp_scalar(p[UU], i, j, k, coeff);
    (*modvar).n_e = rho * Ne_unit + 1e-40;

    Bp[1] = interp_scalar(p[B1], i, j, k, coeff);
    Bp[2] = interp_scalar(p[B2], i, j, k, coeff);
    Bp[3] = interp_scalar(p[B3], i, j, k, coeff);

    V_u[1] = interp_scalar(p[U1], i, j, k, coeff);
    V_u[2] = interp_scalar(p[U2], i, j, k, coeff);
    V_u[3] = interp_scalar(p[U3], i, j, k, coeff);

    VdotV = 0.;
    for (i = 1; i < NDIM; i++)
        for (j = 1; j < NDIM; j++)
            VdotV += g_dd[i][j] * V_u[i] * V_u[j];
    Vfac = sqrt(-1. / g_uu[0][0] * (1. + fabs(VdotV)));
    (*modvar).U_u[0] = -Vfac * g_uu[0][0];
    for (i = 1; i < NDIM; i++)
        (*modvar).U_u[i] = V_u[i] - Vfac * g_uu[0][i];

    lower_index(X, (*modvar).U_u, (*modvar).U_d);

    double Utot = 0;
    for (int i = 0; i < NDIM; i++)
        Utot += (*modvar).U_u[i] * (*modvar).U_d[i];

    UdotBp = 0.;
    for (i = 1; i < NDIM; i++)
        UdotBp += (*modvar).U_d[i] * Bp[i];
    (*modvar).B_u[0] = UdotBp;
    for (i = 1; i < NDIM; i++)
        (*modvar).B_u[i] =
            (Bp[i] + (*modvar).U_u[i] * UdotBp) / (*modvar).U_u[0];

    lower_index(X, (*modvar).B_u, (*modvar).B_d);

    bsq = (*modvar).B_u[0] * (*modvar).B_d[0] +
          (*modvar).B_u[1] * (*modvar).B_d[1] +
          (*modvar).B_u[2] * (*modvar).B_d[2] +
          (*modvar).B_u[3] * (*modvar).B_d[3];

    (*modvar).B = sqrt(bsq) * B_unit + 1e-40;

    beta = uu * (gam - 1.) / 0.5 / bsq;

    beta_trans = 1.;
    b2 = pow(beta / beta_trans, 2);

    trat = 3.;
    two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
    Th_unit = (1.4444444444 - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + trat);

    (*modvar).theta_e =
        (2. / 15.) * (uu / rho) * (PROTON_MASS / ELECTRON_MASS) + 1e-40;

#if (DEBUG)
    if (uu < 0)
        fprintf(stderr, "U %e %e\n", uu, p[UU][i][j][k]);
    ;

    if ((*modvar).theta_e < 0)
        fprintf(stderr, "Te %e\n", (*modvar).theta_e);
    if ((*modvar).B < 0)
        fprintf(stderr, "B %e\n", (*modvar).B);
    if ((*modvar).n_e < 0)
        fprintf(stderr, "Ne %e %e\n", (*modvar).n_e, p[KRHO][i][j][k]);
#endif

    if (bsq / rho > 1. || exp(X[1]) > 50.) {
        (*modvar).n_e = 0;
        return 0;
    }

    return 1;
}

void set_units(double M_unit_) {

    L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    T_unit = L_unit / SPEED_OF_LIGHT;

    RHO_unit = M_unit_ / pow(L_unit, 3);
    U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);

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

    return;
}

void compute_spec_user(struct Camera *intensityfield,
                       double energy_spectrum[num_frequencies][nspec]) {

    return;
}
