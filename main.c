/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *             ___  ___   ___  __________  ___
 *            / _ \/ _ | / _ \/_  __/ __ \/ _ \
 *           / , _/ __ |/ ___/ / / / /_/ / , _/
 *          /_/|_/_/ |_/_/    /_/  \____/_/|_|
 *
 * This program integrates the equations of motion of General Relativity
 * to compute the trajectories of photon bundles (null geodesics); it then
 * performs radiative transfer along these geodesics to compute an image
 * or spectrum. The gravitational field is defined by the metric selected
 * by the user; plasma models can be supplied in the form of GRMHD
 * simulations or analytic models.
 *
 * CONVENTIONS:
 *
 * (Null) geodesics are parametrized by an (affine) parameter called lambda.
 *
 * Metric sign convention: (-,+,+,+)
 *
 * Indices are labeled: "u" (up)   - contravariant index
 *                      "d" (down) - covariant index
 *
 * Examples: U_u[alpha], U_d[alpha], gamma_udd[mu][alpha][beta]
 *
 * A 'ray' (photon bundle position and wave vector) is represented as:
 *
 * photon_u[0] = X_u[0] // Position
 * photon_u[1] = X_u[1]
 * photon_u[2] = X_u[2]
 * photon_u[3] = X_u[3]
 * photon_u[4] = U_u[0] // Wavevector
 * photon_u[5] = U_u[1]
 * photon_u[6] = U_u[2]
 * photon_u[7] = U_u[3]
 *
 * Indices 0, 1, 2, 3 correspond to t, r, theta, phi (Schwarzschild/Kerr).
 */

#include "constants.h"
#include "functions.h"
#include "parameters.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[]) {

    // INPUT FILE
    /////////////

    read_model(argv);

    // INITIALIZE MODEL
    ///////////////////

    // Initialize HARM2D grmhd model
    // Note: this sets the black hole spin 'a'
    init_model();

    // Set constants such as R_ISCO, JANSKY_FACTOR
    // These depend on the black hole spin
    set_constants();

    // INITIALIZE VARIABLES
    ///////////////////////

    // Set the "camera time" (initial coordinate time of rays)
    //    double t_init = 0.; // MUST BE UPDATED
    //    sscanf(argv[1], "%lf", &t_init);

    int steps;

    // Stepsize for constructing the impact parameters alpha, beta
    double stepx = CAM_SIZE_X / (double)IMG_WIDTH;
    double stepy = CAM_SIZE_Y / (double)IMG_HEIGHT;
    double photon_u[8], alpha, beta;
    int x, y, f;

    // MAIN PROGRAM LOOP
    ////////////////////

    // Change this to:
    // initialize_ray(x, y, photon_u);
    // img[y * IMG_WIDTH + x] = integrate_backward(photon_u);

    int num_indices =
        FREQS_PER_DEC * (int)(log10(FREQ_MAX) - log10(FREQ_MIN)) + 1;
    fprintf(stderr, "\nNumber of frequencies to compute: %d\n", num_indices);
    double energy_spectrum[num_indices];
    double frequencies[num_indices];

    double f_x_field[IMG_WIDTH * IMG_HEIGHT];
    double f_y_field[IMG_WIDTH * IMG_HEIGHT];
    double p_field[IMG_WIDTH * IMG_HEIGHT];

    double I_field[IMG_WIDTH * IMG_HEIGHT];
    double Q_field[IMG_WIDTH * IMG_HEIGHT];
    double U_field[IMG_WIDTH * IMG_HEIGHT];
    double V_field[IMG_WIDTH * IMG_HEIGHT];

    for (f = 0; f < num_indices; f++) { // For all frequencies...
        frequencies[f] = FREQ_MIN * pow(10., (double)f / (double)FREQS_PER_DEC);
        energy_spectrum[f] = 0.;
        printf("freq = %+.15e\n", frequencies[f]);
    }

    for (x = 0; x < IMG_WIDTH; x++) { // For all pixel columns...
#pragma omp parallel for default(none) private(f, steps, alpha, beta,          \
                                               photon_u)                       \
    shared(num_indices, energy_spectrum, frequencies, f_x_field, f_y_field,    \
           I_field, Q_field, U_field, V_field, p_field, x, stepx, stepy,       \
           CUTOFF_INNER, IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, CAM_SIZE_Y)        \
        schedule(static, 1)
        for (y = 0; y < IMG_HEIGHT; y++) { // For all pixel rows
            if ((y + x * IMG_HEIGHT) % 100 == 0)
                fprintf(stderr, "current pixel %d of %d\n", y + x * IMG_HEIGHT,
                        IMG_WIDTH * IMG_HEIGHT);
            double *lightpath2 = malloc(9 * max_steps * sizeof(double));

            double *IQUV = malloc(4 * sizeof(double));

            // Compute impact parameters for this pixel
            alpha = -CAM_SIZE_X * 0.5 + (x + 0.5) * stepx;
            beta = -CAM_SIZE_Y * 0.5 + (y + 0.5) * stepy;

            double f_x = 0.;
            double f_y = 0.;
            double p = 0.;

            // INTEGRATE THIS PIXEL'S GEODESIC

            integrate_geodesic(alpha, beta, photon_u, lightpath2, &steps,
                               CUTOFF_INNER);

            // PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
            for (f = 0; f < num_indices; f++) {
                radiative_transfer_polarized(lightpath2, steps,
                                                 frequencies[f], &f_x, &f_y, &p,
                                                 0, IQUV);
                energy_spectrum[f] += IQUV[0];

                /*
                            f_x_field[y * IMG_WIDTH + x] = f_x;
                            f_y_field[y * IMG_WIDTH + x] = f_y;
                            p_field[y * IMG_WIDTH + x] = p;
                */

                I_field[y * IMG_WIDTH + x] = IQUV[0];
                Q_field[y * IMG_WIDTH + x] = IQUV[1];
                U_field[y * IMG_WIDTH + x] = IQUV[2];
                V_field[y * IMG_WIDTH + x] = IQUV[3];
            }
            free(lightpath2);
            free(IQUV);
        }
#pragma omp barrier
    }

    // WRITE OUTPUT FILES
    /////////////////////

    // We open ONE spectrum file and multiple image files (one per frequency)
    //    FILE *spectrum    = fopen("output/spectrum.dat", "w");

    // check if output folder excists
    //if (stat("output", &st) == -1) {
    //    mkdir("output", 0700);
    //}

    for (f = 0; f < num_indices; f++) { // For all frequencies...
        // Create filenames, open files
        char dat_filename[256] = "";
        char vtk_filename[256] = "";
        sprintf(dat_filename, "output/img_data_%e_IQUV.dat", frequencies[f]);
        sprintf(vtk_filename, "output/img_data_%e.vtk", frequencies[f]);
        FILE *imgfile = fopen(dat_filename, "w");
        //    FILE *fp          = fopen(vtk_filename, "w");

        // Write image data to file
        //        write_image_polarized(imgfile, intensityfield[f], f_x_field,
        //        f_y_field, p_field, JANSKY_FACTOR);
        write_image_IQUV(imgfile, I_field, Q_field, U_field, V_field,
                         JANSKY_FACTOR);

        // Close image files
        fclose(imgfile);
        //      fclose(fp);

        // Update spectrum file
        //        fprintf(spectrum, "%+.15e\t%+.15e\n", frequencies[f],
        //        JANSKY_FACTOR * energy_spectrum[f]);
    }

    //  fclose(spectrum);

    // END OF PROGRAM
    /////////////////

    return 0;
}
