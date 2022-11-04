/*
 * model file for HAMR data
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Code written by Wanga Mulaudzi and Jordy Davelaar
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

    fprintf(stderr, "\nStarting read in of HAMR GRMHD data...\n");

    init_grmhd_data(GRMHD_FILE);
}

void readattr(hid_t file_id, const char *attr_name, hid_t mem_type_id,
              void *buf) {
    hid_t attr_id; /* attribute identifier */
    herr_t ret;    /* Return value */
    attr_id = H5Aopen(file_id, attr_name, H5P_DEFAULT);
    ret = H5Aread(attr_id, mem_type_id, buf);
    ret = H5Aclose(attr_id);
}

void readdata(hid_t file_id, const char *attr_name, hid_t mem_type_id,
              hid_t memspace, void *buf) {
    hid_t ds_id;     /* dataset identifier */
    herr_t ret;      /* Return value */
    hid_t dataspace; /* data space identifier */

    ds_id = H5Dopen(file_id, attr_name, H5P_DEFAULT);
    dataspace = H5Dget_space(ds_id); /* dataspace handle */
    ret = H5Dread(ds_id, mem_type_id, memspace, dataspace, H5P_DEFAULT, buf);
    ret = H5Dclose(ds_id);
}

void init_grmhd_data(char *fname) {
    hid_t file_id;    /* File identifier */
    hid_t memspace;   /* memory space identifier */
    hsize_t dimsm[1]; /* memory space dimensions */
    herr_t ret;       /* Return value */
    int RANK_OUT = 1; /* dimension of data array for HDF5 dataset */
    int gridIndex, gridIndex2D;
    double *x1_in, *x2_in, *x3_in, *r_in, *h_in, *ph_in, *RHO_in, *UU_in,
        *U1_in, *U2_in, *U3_in, *B1_in, *B2_in, *B3_in, *gdet_in, *Ucov0_in,
        *Ucon0_in;

    int i, j, z, ieh = 0;
    double x[4], xp[4];
    double rin, hin, phin, gdet, Ucov0, Ucon0;
    double rp, hp, x2temp;
    double dMact = 0, Ladv = 0, t;

    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }

    /* get stanndard HAMR header */
    readattr(file_id, "t", H5T_NATIVE_DOUBLE, &t);
    readattr(file_id, "N1", H5T_NATIVE_INT, &N1);
    readattr(file_id, "N2", H5T_NATIVE_INT, &N2);
    readattr(file_id, "N3", H5T_NATIVE_INT, &N3);
    readattr(file_id, "startx", H5T_NATIVE_DOUBLE, &startx[1]);
    readattr(file_id, "dx", H5T_NATIVE_DOUBLE, &dx[1]);
    readattr(file_id, "a", H5T_NATIVE_DOUBLE, &a);
    readattr(file_id, "gam", H5T_NATIVE_DOUBLE, &gam);
    readattr(file_id, "Rin", H5T_NATIVE_DOUBLE, &Rin);
    readattr(file_id, "Rout", H5T_NATIVE_DOUBLE, &Rout);
    readattr(file_id, "hslope", H5T_NATIVE_DOUBLE, &hslope);
    readattr(file_id, "R0", H5T_NATIVE_DOUBLE, &R0);
    startx[2]=0;
    dx[2]/=2.;

    stopx[0] = 1.;
    stopx[1] = startx[1] + N1 * dx[1];
    stopx[2] = startx[2] + N2 * dx[2];
    stopx[3] = startx[3] + N3 * dx[3];

    double Rh = (1. + sqrt(1. - a * a));


    init_storage();

    /* allocate the memory for dataset */
    x1_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    x2_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    x3_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    r_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    h_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    ph_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    RHO_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    UU_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    // U0_in    = (double *) malloc(N1*N2*N3 * sizeof(double));
    U1_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    U2_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    U3_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    B1_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    B2_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    B3_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    gdet_in = (double *)malloc(N1 * N2 * sizeof(double));
    Ucov0_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));
    Ucon0_in = (double *)malloc(N1 * N2 * N3 * sizeof(double));

    /* memory size of the data */
    dimsm[0] = N1 * N2 * N3;
    memspace = H5Screate_simple(RANK_OUT, dimsm, NULL);

    /* read the datasets */
    readdata(file_id, "x1", H5T_NATIVE_DOUBLE, memspace, &x1_in[0]);
    readdata(file_id, "x2", H5T_NATIVE_DOUBLE, memspace, &x2_in[0]);
    readdata(file_id, "x3", H5T_NATIVE_DOUBLE, memspace, &x3_in[0]);
    readdata(file_id, "r", H5T_NATIVE_DOUBLE, memspace, &r_in[0]);
    readdata(file_id, "h", H5T_NATIVE_DOUBLE, memspace, &h_in[0]);
    readdata(file_id, "ph", H5T_NATIVE_DOUBLE, memspace, &ph_in[0]);
    readdata(file_id, "RHO", H5T_NATIVE_DOUBLE, memspace, &RHO_in[0]);
    readdata(file_id, "UU", H5T_NATIVE_DOUBLE, memspace, &UU_in[0]);
    readdata(file_id, "U1", H5T_NATIVE_DOUBLE, memspace, &U1_in[0]);
    readdata(file_id, "U2", H5T_NATIVE_DOUBLE, memspace, &U2_in[0]);
    readdata(file_id, "U3", H5T_NATIVE_DOUBLE, memspace, &U3_in[0]);
    readdata(file_id, "B1", H5T_NATIVE_DOUBLE, memspace, &B1_in[0]);
    readdata(file_id, "B2", H5T_NATIVE_DOUBLE, memspace, &B2_in[0]);
    readdata(file_id, "B3", H5T_NATIVE_DOUBLE, memspace, &B3_in[0]);
    readdata(file_id, "Ucov0", H5T_NATIVE_DOUBLE, memspace, &Ucov0_in[0]);
    readdata(file_id, "Ucon0", H5T_NATIVE_DOUBLE, memspace, &Ucon0_in[0]);

    fprintf(stderr, "Done!\n");

    /* the memory space of "gdet" is 2D */
    dimsm[0] = N1 * N2;
    memspace = H5Screate_simple(RANK_OUT, dimsm, NULL);

    readdata(file_id, "gdet", H5T_NATIVE_DOUBLE, memspace, &gdet_in[0]);

    /* close HDF5 file */
    ret = H5Fclose(file_id);

    /* find the index for event horizon ridius */
    for (i = 0; i < N1; i++) {
        if (r_in[i * N2 * N3] >= Rh) {
            ieh = i;
            break;
        }
    }

    /* pass the 1D dataset to pointers */
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (z = 0; z < N3; z++) {
                gridIndex = i * N2 * N3 + j * N3 + z;
                gridIndex2D = i * N2 + j;
                x[1] = x1_in[gridIndex];
                x[2] = x2_in[gridIndex];
                x[3] = x3_in[gridIndex];
                rin = r_in[gridIndex];
                hin = h_in[gridIndex];
                phin = ph_in[gridIndex];
                /*
                  H-AMR internal coordinates: x2c = (1+x2)/2
                  --> In grmonty, x2 is treated as x2c
                */
                x2temp = (1.0 + x[2]) / 2.0;
                x[2] = x2temp;
                /* check that we've got the coordinate parameters right */
                coord(i, j, z, xp); // JD: what is this function?
                bl_coord(x, &rp, &hp);
                if (fabs(x[1] - xp[1]) > 1.e5 * x[1] ||
                    fabs(x[2] - xp[2]) > 1.e5 || fabs(x[3] - xp[3]) > 1.e5) {
                    fprintf(stderr, "grid setup error\n");
                    fprintf(stderr,
                            "x[1],xp[1],x[2],xp[2],x[3],xp[3]: %g %g %g %g %g "
                            "%g \n",
                            x[1], xp[1], x[2], xp[2], x[3], xp[3]);
                    fprintf(stderr,
                            "check the internal coordinates, and continue\n");
                    exit(1);
                } else if (fabs(rp - rin) > 1.e-3 * rp ||
                           fabs(hp - hin) > 1.e-5 ||
                           fabs(xp[3] - phin) > 1.e-5) {
                    fprintf(stderr, "grid setup error\n");
                    fprintf(
                        stderr,
                        "rp,r, (rp-r)/r, hp,h,php,ph: %g %g %g %g %g %g %g \n",
                        rp, rin, (rp - rin) / rp, hp, hin, xp[3], phin);
                    fprintf(stderr, "edit R0, hslope, compile, and continue\n");
                    exit(1);
                }

                /* Since x2c = (1+x2)/2, the vectors in x2 direction should
                   be corrected. Or, we need to correct the theta correction
                   term in gcov_func  (pi -> pi/2)
                */
                p[KRHO][i][j][z] = RHO_in[gridIndex];
                p[UU][i][j][z] = UU_in[gridIndex];
                p[U1][i][j][z] = U1_in[gridIndex];
                p[U2][i][j][z] = U2_in[gridIndex] / 2.;
                p[U3][i][j][z] = U3_in[gridIndex];
                p[B1][i][j][z] = B1_in[gridIndex];
                p[B2][i][j][z] = B2_in[gridIndex] / 2.;
                p[B3][i][j][z] = B3_in[gridIndex];
                gdet = gdet_in[gridIndex2D];
                Ucov0 = Ucov0_in[gridIndex];
                Ucon0 = Ucon0_in[gridIndex];

                /* check accretion rate */
                if (i == ieh)
                    dMact += gdet * p[KRHO][i][j][z] * p[U1][i][j][z] * Ucon0;
                if (i >= 20 && i < 40)
                    Ladv +=
                        gdet * p[UU][i][j][z] * p[U1][i][j][z] * Ucon0 * Ucov0;
            }
        }
    }

    /* deallocate memories */
    free(x1_in);
    free(x2_in);
    free(x3_in);
    free(r_in);
    free(h_in);
    free(ph_in);
    free(RHO_in);
    free(UU_in);
    // free(U0_in);
    free(U1_in);
    free(U2_in);
    free(U3_in);
    free(B1_in);
    free(B2_in);
    free(B3_in);
    free(gdet_in);
    free(Ucov0_in);
    free(Ucon0_in);

    // bias_norm /= V;
    // dMact *= dx[3] * dx[2];
    /* since dx[2] was rearranged by dx[2] = dx[2]/2 while using gdet from the
      data, the accretion rate should be multiplied by 2 */
    dMact *= dx[3] * dx[2] * 2;
    // dMact /= 21.;
    Ladv *= dx[3] * dx[2];
    Ladv /= 21.;
}

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th) {

    *r = exp(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

    return;
}

void coord(int i, int j, int k, double *X) {

    /* returns zone-centered values for coordinates */
    X[0] = startx[0];
    X[1] = startx[1] + (i + 0.5) * dx[1];
    X[2] = startx[2] + (j + 0.5) * dx[2];
    X[3] = startx[3] + (k) * dx[3];

    return;
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
