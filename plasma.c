/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka
 *
 */

#include <math.h>
#include "functions.h"
#include "constants.h"
#include "parameters.h"
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

// Return electron density at position X_u
void disk_velocity(double X_u[4], double Uplasma_u[4]){
    double omega, r2, r32, rho2, cth2, sth2;
    double g_tt, g_tp, g_pp;
    double Ucon_bl[DIM];
    double m1oA2, AA;

    double spot_ne, spot_ne0, Rspot;
    double xspot[DIM];//curent position of a spot center in KS'
    double th_spot, r_spot;
    double xx;//, ux;
    double P;// X_cov[DIM], xspot_cov[DIM];
    double xc, yc, zc;
    double xs, ys, zs;

    double r  = logscale ? exp(X_u[1]) : X_u[1];
    double th = X_u[2];

    //backgound density
    rho2 = r * r;

    //Keplerian velocity of plasma
    sth2  = sin(th) * sin(th);
    cth2  = cos(th) * cos(th);
    r2    = r * r;
    r32   = pow(r, 1.5);
    rho2  = r2 + a * a * cth2;
    g_tt  = -(1. - 2. * r / rho2);
    g_tp  = -2. * a * r * sth2 / rho2;
    g_pp  = (r2 + a * a + 2 * a * a * r * sth2 / rho2) * sth2;
    omega = 1. / (r32 + a);
    m1oA2 = g_tt + 2. * omega * g_tp + omega * omega * g_pp;
    AA    = sqrt(-1. / m1oA2);

    // Formula for Keplerian velocity in BL metric, notice that has to be transformedinto MBL,Ks or MKS
    Ucon_bl[0] = AA;
    Ucon_bl[1] = 0.;
    Ucon_bl[2] = 0.;
    Ucon_bl[3] = AA * omega;

    int i;

    // In case of KS/MKS coords, convert this four-vector into corresponding coordinate system
    #if(metric == KS || metric == MKS)
        double BLcoords[8], KScoords[8];
        LOOP_i{
            BLcoords[i] = X_u[i];
            BLcoords[i+4] = Ucon_bl[i];
        }
        BL_to_KS_u(BLcoords, KScoords);
        LOOP_i Ucon_bl[i] = KScoords[i+4];
    #endif

    // Put plasma velocity into u_u[4]
    LOOP_i Uplasma_u[i] = Ucon_bl[i];
}

// Return the magnetic field intensity at position X_u
// B strength defined through beta plasma parameter
double B_fun(double X_u[4], double ne){
    double r    = logscale ? exp(X_u[1]) : X_u[1];
    double beta = 10.; //plasma parameter
    return sqrt(8. * M_PI / beta * ne * MPCL2 * 2. / 12. / r);
}

// Return the value for Te at position X_u
double THETA_e(double X_u[4]){
    double r = logscale ? exp(X_u[1]) : X_u[1];
    return THETA_e_0 * pow(r, -0.84);
}
