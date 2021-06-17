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


#include "raptor_harm_model.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "functions.h"

#include "parameters.h"

void init_model()
{
    /* find dimensional quantities from black hole
     mass and its accretion rate */
    set_units(M_UNIT);

    fprintf(stderr, "getting simulation data...\n");

    init_harm_data(GRMHD_FILE);
}


void init_harm_data(char *fname)
{
    FILE *fp;
    double x[4];
//    double rp, hp, V, two_temp_gam; // These appear to be unnecessary
//    double dV;
    int i, j, k;

    /* header variables not used except locally */
    double t, tf, cour, DTd, DTl, DTi, dt;
    int nstep, DTr, dump_cnt, image_cnt, rdump_cnt, lim, failed;
    double r, h, divb, vmin, vmax, gdet;
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    fprintf(stderr, "%s\n", fname);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        fprintf(stderr, "can't open sim data file\n");
        exit(1);
    } else {
        fprintf(stderr, "successfully opened %s\n", fname);
    }

    /* get standard HARM header */
    fscanf(fp, "%lf ", &t);
    fscanf(fp, "%d ", &N1);
    fscanf(fp, "%d ", &N2);
    fscanf(fp, "%lf ", &startx[1]);
    fscanf(fp, "%lf ", &startx[2]);
    fscanf(fp, "%lf ", &dx[1]);
    fscanf(fp, "%lf ", &dx[2]);
    fscanf(fp, "%lf ", &tf);
    fscanf(fp, "%d ", &nstep);
    fscanf(fp, "%lf ", &a);
    fscanf(fp, "%lf ", &gam);
    fscanf(fp, "%lf ", &cour);
    fscanf(fp, "%lf ", &DTd);
    fscanf(fp, "%lf ", &DTl);
    fscanf(fp, "%lf ", &DTi);
    fscanf(fp, "%d ", &DTr);
    fscanf(fp, "%d ", &dump_cnt);
    fscanf(fp, "%d ", &image_cnt);
    fscanf(fp, "%d ", &rdump_cnt);
    fscanf(fp, "%lf ", &dt);
    fscanf(fp, "%d ", &lim);
    fscanf(fp, "%d ", &failed);
    fscanf(fp, "%lf ", &Rin);
    fscanf(fp, "%lf ", &Rout);
    fscanf(fp, "%lf ", &hslope);
    fscanf(fp, "%lf ", &R0);
//printf("hslope = %lf\n", hslope);

    // nominal non-zero values for axisymmetric simulations
    startx[0] = 0.;
    startx[2]=startx[2];// *M_PI;
    startx[3] = 0.;
    dx[2]=dx[2];// *M_PI;
    stopx[0] = 1.;
    stopx[1] = startx[1] + N1 * dx[1];
    stopx[2] = startx[2] + N2 * dx[2];
    stopx[3] = 2. * M_PI;

    fprintf(stderr, "Sim range x1, x2:  %g %g, %g %g\n", startx[1],
            stopx[1], startx[2], stopx[2]);

    dx[0] = 1.;
    dx[3] = 2. * M_PI;

    // Allocate storage for all model size dependent variables
    init_storage();

//    two_temp_gam =
 //   0.5 * ((1. + 2. / 3. * (TP_OVER_TE + 1.) / (TP_OVER_TE + 2.)) +
  //         gam);
    Thetae_unit = (gam - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + TP_OVER_TE);

    dMact = 0.;
    Ladv = 0.;
//    V = 0.;
//    dV = dx[1] * dx[2] * dx[3];
    for (k = 0; k < N1 * N2; k++) {
        j = k % N2;
        i = (k - j) / N2;
        fscanf(fp, "%lf %lf %lf %lf", &x[1], &x[2], &r, &h);

        // Check that we've got the coordinate parameters right
  /*      bl_coord(x, &rp, &hp);
        if (fabs(rp - r) > 1.e-5 * rp || fabs(hp - h) > 1.e-5) {
            fprintf(stderr, "grid setup error\n");
            fprintf(stderr, "rp,r,hp,h: %g %g %g %g\n",
                    rp, r, hp, h);
            fprintf(stderr,
                    "edit R0, hslope, compile, and continue\n");
            exit(1);
        }
*/
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf",
               &p[KRHO][i][j],
               &p[UU][i][j],
               &p[U1][i][j],
               &p[U2][i][j],
               &p[U3][i][j],
               &p[B1][i][j], &p[B2][i][j], &p[B3][i][j]);
//	fprintf(stderr, "%lf\n", p[KRHO][i][j]);

        fscanf(fp, "%lf", &divb);

        fscanf(fp, "%lf %lf %lf %lf",
               &Ucon[0], &Ucon[1], &Ucon[2], &Ucon[3]);
        fscanf(fp, "%lf %lf %lf %lf", &Ucov[0],
               &Ucov[1], &Ucov[2], &Ucov[3]);
        fscanf(fp, "%lf %lf %lf %lf", &Bcon[0],
               &Bcon[1], &Bcon[2], &Bcon[3]);
        fscanf(fp, "%lf %lf %lf %lf", &Bcov[0],
               &Bcov[1], &Bcov[2], &Bcov[3]);
        fscanf(fp, "%lf ", &vmin);
        fscanf(fp, "%lf ", &vmax);
        fscanf(fp, "%lf ", &vmin);
        fscanf(fp, "%lf ", &vmax);
        fscanf(fp, "%lf\n", &gdet);

        // Check accretion rate
        if (i <= 20)
            dMact += gdet * p[KRHO][i][j] * Ucon[1];
        if (i >= 20 && i < 40)
            Ladv += gdet * p[UU][i][j] * Ucon[1] * Ucov[0];
    }

    dMact *= dx[3] * dx[2];
    dMact /= 21.;
    Ladv *= dx[3] * dx[2];
    Ladv /= 21.;
    fprintf(stderr, "dMact: %g, Ladv: %g\n", dMact, Ladv);
    fprintf(stderr, "Done reading data\n");
}


// TRANSFORMATION FUNCTIONS
///////////////////////////

// WARNING: these are not yet listed in functions.h and are only meant for use
// by other functions in this file.

// Returns the value of f(Xg2) given some value for Xr2. For the correct Xg2,
// we have f(Xg2) = 0.
double f_Xg2(double Xg2, double Xr2){
    return M_PI * Xg2 + 0.5 * (1. - hslope) * sin(2. * M_PI * Xg2) - Xr2;
}

// Returns the value of f'(Xg2).
double f_primed_Xg2(double Xg2){
    return M_PI + M_PI * (1. - hslope) * cos(2. * M_PI * Xg2);
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Xg2_1 based on the input estimate Xg2_0.
double NR_stepX(double Xg2_0, double Xr2){
    double fprime = f_primed_Xg2(Xg2_0);

    if(fabs(fprime) < 1.e-9)
        printf("fprime = %+.15e\n", fprime);

    return Xg2_0 - f_Xg2(Xg2_0, Xr2) / f_primed_Xg2(Xg2_0);
}

// Returns the value of f(Ug2) given some value for Ur2. For the correct Ug2,
// we have f(Ug2) = 0.
double f_Ug2(double Ug2, double Ur2, double Xg2){
    return M_PI * Ug2 * (1. + (1. - hslope) * cos(2. * M_PI * Xg2)) - Ur2;
}

// Returns the value of f'(Ug2).
double f_primed_Ug2(double Ug2, double Xg2){
    return M_PI * (1. + (1. - hslope) * cos(2. * M_PI * Xg2));
}

// This function does "one Newton-Raphson step", i.e. it returns the NEW,
// "better" estimate Ug2_1 based on the input estimate Ug2_0.
double NR_stepU(double Ug2_0, double Ur2, double Xg2){
    double fprime = f_primed_Ug2(Ug2_0, Xg2);

    if(fabs(fprime) < 1.e-9)
        printf("fprime = %+.15e\n", fprime);

    return Ug2_0 - f_Ug2(Ug2_0, Ur2, Xg2) / f_primed_Ug2(Ug2_0, Xg2);
}

// Given the X2 coordinate in RAPTOR's convention, Xr2, we compute and return
// an estimate for the corresponding coordinate in HARM2D's convention, Xg2.
double Xg2_approx_rand(double Xr2){
    double Xg2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
    double Xg2_prev = 1.e-15;     // Keeps track of previous estimate to converge
    double tolerance = 1.e-9; // Maximum error
    int steps = 0;
    int maxsteps = 100;

    int count = 0;

    // Main loop
    while (fabs(Xg2_current - Xg2_prev) > tolerance){
        Xg2_current = (double) rand() / (double)RAND_MAX;
        //Xg2_current = 1.e-16;
        steps = 0;
        count++;

        while(steps < maxsteps && fabs(Xg2_current - Xg2_prev) > tolerance){
            Xg2_prev = Xg2_current;
            Xg2_current = NR_stepX(Xg2_current, Xr2);
            steps++;
        }
    }

    // Clamp output value between 0 and 1
    return fmin(1., fmax(Xg2_current, 0.));
}

// Given the U2 coordinate in RAPTOR's convention, Ur2, we compute and return
// an estimate for the corresponding vector component in HARM2D's convention, Ug2.
double Ug2_approx_rand(double Ur2, double Xg2){
    double Ug2_current = 0.1; // Initial guess; reasonable b/c Xg2 E [0, 1]
    double Ug2_prev = 1.e-15;     // Keeps track of previous estimate to converge
    double tolerance = 1.e-9; // Maximum error
    int steps = 0;
    int maxsteps = 100;

    int count = 0;

    // Main loop
    while (fabs(Ug2_current - Ug2_prev) > tolerance){
        Ug2_current = (double) rand() / (double)RAND_MAX;
        steps = 0;
        count++;

        while(steps < maxsteps && fabs(Ug2_current - Ug2_prev) > tolerance){
            Ug2_prev = Ug2_current;
            Ug2_current = NR_stepU(Ug2_current, Ur2, Xg2);
            steps++;
        }
    }

    return Ug2_current;
}


// Current metric: modified Kerr-Schild, squashed in theta
// to give higher resolution at the equator


#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)

/* mnemonics for dimensional indices */
#define TT      0
#define RR      1
#define TH      2
#define PH      3

/* return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]) + R0;
	*th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);

	return;
}

void gcon_func(double *X, double gcon[][NDIM])
{

	int k, l;
	double sth, cth, irho2;
	double r, th;
	double hfac;
	/* required by broken math.h */
	//void sincos(double in, double *sth, double *cth);

	DLOOP gcon[k][l] = 0.;

	bl_coord(X, &r, &th);

	sth=sin(th);//, &sth, &cth);
	sth = fabs(sth) + 1.e-9;

	irho2 = 1. / (r * r + a * a * cth * cth);

	// transformation for Kerr-Schild -> modified Kerr-Schild 
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);

	gcon[TT][TT] = -1. - 2. * r * irho2;
	gcon[TT][1] = 2. * irho2;

	gcon[1][TT] = gcon[TT][1];
	gcon[1][1] = irho2 * (r * (r - 2.) + a * a) / (r * r);
	gcon[1][3] = a * irho2 / r;

	gcon[2][2] = irho2 / (hfac * hfac);

	gcon[3][1] = gcon[1][3];
	gcon[3][3] = irho2 / (sth * sth);
}


void gcov_func(double *X, double gcov[][NDIM])
{
	int k, l;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;
	/* required by broken math.h */
	void sincos(double th, double *sth, double *cth);

	DLOOP gcov[k][l] = 0.;

	bl_coord(X, &r, &th);

	sth=sin(th);//, &sth, &cth);
	sth = fabs(sth) + 1.e-9;
	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	/* transformation for Kerr-Schild -> modified Kerr-Schild */
	tfac = 1.;
	rfac = r - R0;
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}

#undef TT
#undef RR
#undef TH
#undef PH


void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

	ucov[0] = Gcov[0][0] * ucon[0]
	    + Gcov[0][1] * ucon[1]
	    + Gcov[0][2] * ucon[2]
	    + Gcov[0][3] * ucon[3];
	ucov[1] = Gcov[1][0] * ucon[0]
	    + Gcov[1][1] * ucon[1]
	    + Gcov[1][2] * ucon[2]
	    + Gcov[1][3] * ucon[3];
	ucov[2] = Gcov[2][0] * ucon[0]
	    + Gcov[2][1] * ucon[1]
	    + Gcov[2][2] * ucon[2]
	    + Gcov[2][3] * ucon[3];
	ucov[3] = Gcov[3][0] * ucon[0]
	    + Gcov[3][1] * ucon[1]
	    + Gcov[3][2] * ucon[2]
	    + Gcov[3][3] * ucon[3];

	return;
}

// Get the fluid parameters - IN THE PLASMA FRAME?
void get_fluid_params(double X[NDIM], double *Ne,
                      double *Thetae, double *B, double Bcon[NDIM], double Ucon[NDIM], int *IN_VOLUME)
{
    int i, j;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;
    double gcon[NDIM][NDIM], gcov[NDIM][NDIM], coeff[4];
    double interp_scalar(double **var, int i, int j, double del[4]);
    double Ucov[NDIM], Bcov[NDIM];


#if(metric == MKS)
    double Xtemp = X[2]; // Remember the old theta - we want to go back to MKS!
    // TRANSFORM FROM LKS TO MKS (as defined by Gammie et al) COORDINATES
    X[2] = Xg2_approx_rand(X[2]); // We only transform theta - r is already exponential and R0 = 0
#endif

//    metric_uu(X, gcon);
//    metric_dd(X, gcov);
gcon_func(X, gcon);
gcov_func(X, gcov); // Gammies metric

    *IN_VOLUME = 1;

    if (X[1] < startx[1] ||
        X[1] > stopx[1] || X[2] < startx[2] || X[2] > stopx[2]) {
        *Ne = 0.;

#if(metric == MKS)
        X[2] = Xtemp; // Need to go back to LKS!
#endif

        *IN_VOLUME = 0;
        return;
    }

double smalll = 1.e-20;

    // conversie van coordinaat naar cel i,j
    Xtoij(X, &i, &j, del);

    //coef nodig voor je interpolatie
    coeff[0] = (1. - del[1]) * (1. - del[2]);
    coeff[1] = (1. - del[1]) * del[2];
    coeff[2] = del[1] * (1. - del[2]);
    coeff[3] = del[1] * del[2];
    //inteprolatie van je primitieve variabelen
    rho = interp_scalar(p[KRHO], i, j, coeff) + smalll;
    uu = interp_scalar(p[UU], i, j, coeff);

    //bepalen van de plasma number density en electron temperatuur
    *Ne = rho * Ne_unit;
    *Thetae = uu / rho * Thetae_unit + smalll;

    Bp[1] = interp_scalar(p[B1], i, j, coeff);
    Bp[2] = interp_scalar(p[B2], i, j, coeff);
    Bp[3] = interp_scalar(p[B3], i, j, coeff);

    Vcon[1] = interp_scalar(p[U1], i, j, coeff);
    Vcon[2] = interp_scalar(p[U2], i, j, coeff);
    Vcon[3] = interp_scalar(p[U3], i, j, coeff);

    //reconstrueren van de 4 vectoren
    // Get Ucov
    VdotV = 0.;
    for (i = 1; i < NDIM; i++)
        for (j = 1; j < NDIM; j++)
            VdotV += gcov[i][j] * Vcon[i] * Vcon[j];
    Vfac = sqrt(-1. / gcon[0][0] * (1. + fabs(VdotV)));
    Ucon[0] = -Vfac * gcon[0][0];
    for (i = 1; i < NDIM; i++)
        Ucon[i] = Vcon[i] - Vfac * gcon[0][i];

//    lower_index(X, Ucon, Ucov);
    lower(Ucon, gcov, Ucov); // Gammie's lowering function

    // Get B and Bcov
    UdotBp = 0.;
    for (i = 1; i < NDIM; i++)
        UdotBp += Ucov[i] * Bp[i];
    Bcon[0] = UdotBp;
    for (i = 1; i < NDIM; i++)
        Bcon[i] = (Bp[i] + Ucon[i] * UdotBp) / Ucon[0];

//    lower_index(X, Bcon, Bcov);
    lower(Bcon, gcov, Bcov); // Gammie's lowering function

    //sterkte van het magneetveld
    *B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
              Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]) * B_unit + smalll;

    double Bsq = (Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] +
              Bcon[2] * Bcov[2] + Bcon[3] * Bcov[3]);

    char diskstring[100] = "DISK";
    strcpy(diskstring, "DISK");
/*
    if (strcmp(TEMP_MODEL, diskstring) == 0){
        double beta_trans=1.;

        double b2=pow(uu*(4./3.-1.) / (0.5 *(Bsq +1e-40) *beta_trans),2);

        double Rhigh=1.;
        double Rlow=1.;
        double trat = Rhigh * b2/(1. + b2) + Rlow /(1. + b2);

        double two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat +2.)) + 4./3.);

        Thetae_unit = (two_temp_gam - 1.) * (MPoME) / (1. + trat);

        *Thetae = (uu/(rho+1e-20))*Thetae_unit;
    }
*/
//    if(Bsq / p[KRHO][i][j] > 1.){
//	*Thetae=30.;
//   }

    // MICHAEL JANSSEN TEMPERATURE MODEL
//    double dummy_1 = 2. * (gam - 1.) * p[UU][i][j] / (Bsq + 1.e-20);
//    double dummy_2 = dummy_1 * dummy_1;
//    double dummy_3 = 1. + dummy_2;
//    double dummy_4 = 100. * dummy_2 / dummy_3  + 1. / dummy_3;
//    double dummy_5 = 0.5 * ((1. + 2. / 3. * (dummy_4 + 1.) / (dummy_4 + 2.)) + gam);
//    if (TEMP_MODEL == 2){
//        *Thetae = p[UU][i][j] / (*Ne) * Ne_unit * (dummy_5 - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + dummy_4);
//        if(Bsq / p[KRHO][i][j] > 1.){
//            *Thetae = 0.;
//        }
  //  }

    if (SPHERICAL_ACC){
        double r, theta;
        bl_coord(X, &r, &theta);
        *B = sqrt(8. * M_PI * rho / r) * B_unit; // MODIFIED TO ADD B FIELD TO BONDI ACCRETION
    }

//printf("\nIN RAPTORMODEL u dot b = %+.15e", Ucon[0] * Bcov[0] + Ucon[1] * Bcov[1] + Ucon[2] * Bcov[2] + Ucon[3] * Bcov[3]);
//if(fabs(Ucon[0] * Ucov[0] + Ucon[1] * Ucov[1] + Ucon[2] * Ucov[2] + Ucon[3] * Ucov[3] + 1.) > 1.e-9)
//printf("\nIN RAPTORMODEL u dot u = %+.15e", Ucon[0] * Ucov[0] + Ucon[1] * Ucov[1] + Ucon[2] * Ucov[2] + Ucon[3] * Ucov[3]);

#if(metric == MKS)
// Transform from MKS back to LKS for further use in RAPTOR
Bcon[2] = M_PI * Bcon[2] * (1. + (1. - hslope) * cos(2. * M_PI * X[2]));
Ucon[2] = M_PI * Ucon[2] * (1. + (1. - hslope) * cos(2. * M_PI * X[2]));
X[2] = Xtemp;
#endif

if(0){
printf("\nBcon[0] = %+.15e\n", Bcon[0]);
printf("Bcon[1] = %+.15e\n", Bcon[1]);
printf("Bcon[2] = %+.15e\n", Bcon[2]);
printf("Bcon[3] = %+.15e\n", Bcon[3]);

printf("Ucon[0] = %+.15e\n", Ucon[0]);
printf("Ucon[1] = %+.15e\n", Ucon[1]);
printf("Ucon[2] = %+.15e\n", Ucon[2]);
printf("Ucon[3] = %+.15e\n", Ucon[3]);

printf("\nBcov[0] = %+.15e\n", Bcov[0]);
printf("Bcov[1] = %+.15e\n", Bcov[1]);
printf("Bcov[2] = %+.15e\n", Bcov[2]);
printf("Bcov[3] = %+.15e\n", Bcov[3]);

printf("Ucov[0] = %+.15e\n", Ucov[0]);
printf("Ucov[1] = %+.15e\n", Ucov[1]);
printf("Ucov[2] = %+.15e\n", Ucov[2]);
printf("Ucov[3] = %+.15e\n", Ucov[3]);

printf("B = %+.15e\n", *B);
printf("Ne = %+.15e\n", *Ne);
printf("Thetae = %+.15e\n", *Thetae);
}


}

void set_units(double M_unit_)
{
//	double MBH;

	/* set black hole mass */
	/** could be read in from file here,
	    along with M_unit and other parameters **/
//	MBH = 4.e6;

	/** input parameters appropriate to Sgr A* **/
	//double BH_MASS = MBH * MSUN;

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

double interp_scalar(double **var, int i, int j, double coeff[4])
{

    double interp;

    interp =
    var[i][j] * coeff[0] +
    var[i][j + 1] * coeff[1] +
    var[i + 1][j] * coeff[2] + var[i + 1][j + 1] * coeff[3];

    return interp;
}

void Xtoij(double X[NDIM], int *i, int *j, double del[NDIM])
{

    *i = (int) ((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
    *j = (int) ((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;

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

    return;
}

// ALLOCATION STUFF BELOW HERE
//////////////////////////////

static void *malloc_rank1(int n1, int alloc_size)
{
    void *A;

    if ((A = (void *) malloc(n1 * alloc_size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

/*
static void **malloc_rank2(int n1, int n2, int alloc_size)
{

    void **A;
    int i;

    if ((A = (void **) malloc(n1 * sizeof(void *))) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank2\n");
        exit(124);
    }

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, alloc_size);
    }

    return A;
}
*/

static double **malloc_rank2_cont(int n1, int n2)
{

    double **A;
    double *space;
    int i;

    space = malloc_rank1(n1 * n2, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++)
        A[i] = &(space[i * n2]);

    return A;
}

void init_storage()
{
    int i;

    p = malloc_rank1(NPRIM, sizeof(double *));
    for (i = 0; i < NPRIM; i++)
        p[i] = (double **) malloc_rank2_cont(N1, N2);

    return;
}
