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

double radiative_transfer_unpolarized(double *lightpath, int steps,
                                      double frequency, double *IQUV) {

    int path_counter;
    double pitch_ang, nu_p;

    double X_u[4], k_d[4], k_u[4], dl_current;
    double jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV;

    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm
    double g_dd[4][4], g_uu[4][4];

    double Icurrent = IQUV[0];

    double dtau_old = 0;

    struct GRMHD modvar;
    modvar.B = 0;
    modvar.n_e = 0.;
    modvar.theta_e = 0;

    LOOP_i {
        modvar.B_u[i] = 0;
        modvar.U_u[i] = 0;
        modvar.B_d[i] = 0;
        modvar.U_d[i] = 0;
    }
    modvar.igrid_c = -1;

    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        metric_dd(X_u, g_dd);
        metric_uu(X_u, g_uu);

        if (get_fluid_params(X_u, &modvar)) {
            // Obtain pitch angle: still no units (geometric)
            lower_index(X_u, k_u, k_d);
            pitch_ang = pitch_angle(X_u, k_u, modvar.B_u, modvar.U_u);

            dl_current *= (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT) /
                          (PLANCK_CONSTANT * frequency);

            // CGS UNITS USED FROM HERE ON OUT
            //////////////////////////////////
            //  if(tau[f] < log(1000.) ) {

            // Scale the wave vector to correct energy
            LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                             (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            // lower the index of the wavevector
            lower_index(X_u, k_u, k_d);

            // Compute the photon frequency in the plasma frame:
            nu_p = freq_in_plasma_frame(modvar.U_u, k_d);
            // Obtain emission coefficient in current plasma conditions
            evaluate_coeffs_single(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI, &aQ, &aU,
                            &aV, nu_p, modvar.theta_e, modvar.n_e, modvar.B,
                            pitch_ang);

            double C = Rg * PLANCK_CONSTANT /
                       (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            double dtau = (aI * dl_current * C + dtau_old);
            double K_inv = aI;
            double j_inv = jI;

            if (jI == jI) { // I_current += exp(-tau) * j_nu /
                            // nu_p / nu_p * dl_current * C;
                //         I_current += dI; // Old way of integrating
                double Ii = Icurrent;
                double S = j_inv / K_inv;
                if (K_inv == 0)
                    Icurrent = Ii;
                else if (dtau < 1.e-5)
                    Icurrent = Ii - (Ii - S) * (0.166666667 * dtau *
                                                (6. - dtau * (3. - dtau)));
                else {
                    double efac = exp(-dtau);
                    Icurrent = Ii * efac + S * (1. - efac);
                }
            }
            dtau_old = 0;

            //      }
        }

        /*
           if(0 && (j_nu != j_nu || isnan(j_nu))){
           printf("NaN emissivity! X_u[2], theta = %+.15e %+.15e\n",
           X_u[2],M_PI*X_u[2] + 2*0.35/M_PI *
           sin(2*M_PI*X_u[2])*atan2(2.*(log(50.)-X_u[1]),1.)); printf("NaN
           emissivity! expX_u[1], innercutoff = %+.15e %e\n",
           exp(X_u[1]),cutoff_inner); printf("NaN emissivity! pitch_angle =
           %+.15e\n", pitch_ang); printf("NaN emissivity! B = %+.15e\n", B);
           printf("NaN emissivity! THETAe = %+.15e\n", THETA_e);
           printf("NaN emissivity! nu_plasmaframe = %+.15e\n", nu_p);
           printf("NaN emissivity! n_e = %+.15e\n", n_e);
           printf("NaN emissivity! j_nu = %+.15e\n", j_nu);
           printf("NaN emissivity! tau = %+.15e\n", tau);
           printf("NaN emissivity! U dot U = %+.15e\n\n", inner_product(X_u,
           Uplasma_u, Uplasma_u)); printf("NaN emissivity! k dot k =
           %+.15e\n\n", inner_product(X_u, k_u, k_u));
           }
         */
    }

    IQUV[0] = Icurrent * pow(frequency, 3.);

    return 1;
}
