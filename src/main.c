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

#include "functions.h"
#include "parameters.h"

int main(int argc, char *argv[]) {

    // INPUT FILE
    /////////////

    fprintf(stderr, "=====================================================\n");
    fprintf(stderr, "             ___  ___   ___  __________  ___         \n");
    fprintf(stderr, "            / _ \/ _ | / _ \/_  __/ __ \/ _ \        \n");
    fprintf(stderr, "           / , _/ __ |/ ___/ / / / /_/ / , _/        \n");
    fprintf(stderr, "          /_/|_/_/ |_/_/    /_/  \____/_/|_|         \n");
    fprintf(stderr, "=====================================================\n");

    read_model(argv);

    // INITIALIZE MODEL
    ///////////////////

    // Initialize HARM2D grmhd model
    // Note: this sets the black hole spin 'a'
    init_model();

    // Set constants such as R_ISCO, JANSKY_FACTOR
    // These depend on the black hole spin
    set_constants();

    // MAIN PROGRAM LOOP
    ////////////////////

    // Change this to:
    // initialize_ray(x, y, photon_u);
    // img[y * IMG_WIDTH + x] = integrate_backward(photon_u);

    fprintf(stderr, "\nNumber of frequencies to compute: %d\n",
            num_frequencies);
    double energy_spectrum[num_frequencies];
    double frequencies[num_frequencies];

    struct Camera *intensityfield;

    init_camera(&intensityfield);

    for (int f = 0; f < num_frequencies; f++) { // For all frequencies...
        frequencies[f] = FREQ_MIN * pow(10., (double)f / (double)FREQS_PER_DEC);
        energy_spectrum[f] = 0.;
        printf("freq = %+.15e\n", frequencies[f]);
    }

    int block = 0;

    while (block < tot_blocks) {
        if (block % (10) == 0)
            fprintf(stderr, "block %d of total %d\n", block, tot_blocks);

        calculate_image_block(&intensityfield[block], energy_spectrum,
                              frequencies);

        if (refine_block(intensityfield[block])) {
            add_block(&intensityfield, block);
        } else {
            block++;
        }
    }

    compute_spec(intensityfield, energy_spectrum);

    // WRITE OUTPUT FILES
    /////////////////////

    output_files(intensityfield, energy_spectrum, frequencies);

    fprintf(stderr, "\nFinished writing files.\n");

    // FREE ALLOCATED POINTERS
    //////////////////////////

    free(intensityfield);

    // END OF PROGRAM
    /////////////////

    return 0;
}
