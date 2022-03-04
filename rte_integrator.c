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

double radiative_transfer(double *lightpath, int steps, double frequency) {
    int path_counter;
    double I_current = 0.;
    double j_nu = 0.;
    double pitch_ang, nu_p, dl_current;
    int i;
    double X_u[4], k_u[4], k_d[4];
    double Rg = GGRAV * MBH / SPEED_OF_LIGHT / SPEED_OF_LIGHT; // Rg in cm

    double a_nu = 0.;
    double dtau_old = 0;

    struct GRMHD modvar;

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

        // Check whether the ray is currently in the GRMHD simulation volume
        if (get_fluid_params(X_u, modvar)) {
            // Obtain pitch angle: still no units (geometric)
            pitch_ang = pitch_angle(X_u, k_u, modvar.B_u, modvar.U_u);

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
            nu_p = freq_in_plasma_frame(modvar.U_u, k_d);

            // Obtain emission coefficient in current plasma conditions
            j_nu = emission_coeff_THSYNCHAV(modvar.B, modvar.theta_e, nu_p,
                                            modvar.n_e);

            // Obtain absorption coefficient
            if (ABSORPTION) {
                a_nu = absorption_coeff_TH(j_nu, nu_p, modvar.theta_e);
            }

            // Constant used in integration (to produce correct units)
            double C = Rg * PLANCK_CONSTANT /
                       (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

            double dtau = (nu_p * a_nu * dl_current * C + dtau_old);
            double K_inv = (nu_p * a_nu);
            double j_inv = (j_nu / (nu_p * nu_p));

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
