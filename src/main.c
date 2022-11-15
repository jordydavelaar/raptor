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

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

int main(int argc, char *argv[]) {

    // INPUT FILE
    /////////////

    fprintf(stderr,
            "__________    _____ _____________________________ __________\n");
    fprintf(stderr, "\\______   \\  /  _  \\\\______   \\__    ___/\\_____  "
                    "\\\\______   \\\n");
    fprintf(
        stderr,
        " |       _/ /  /_\\  \\|     ___/ |    |    /   |   \\|       _/\n");
    fprintf(
        stderr,
        " |    |   \\/    |    \\    |     |    |   /    |    \\    |   \\\n");
    fprintf(
        stderr,
        " |____|_  /\\____|__  /____|     |____|   \\_______  /____|_  /\n");
    fprintf(
        stderr,
        "        \\/         \\/                            \\/       \\/  \n");

    fprintf(stderr, "\nRunning RAPTOR v1.0 in");

    if (POL)
        fprintf(stderr, " polarized mode!\n");
    else
        fprintf(stderr, " unpolarized mode!\n");

    fprintf(stderr, "\nInitializing...\n");
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

    fprintf(stderr, "\nNumber of frequencies to compute: %d\n",
            num_frequencies);
    double energy_spectrum[num_frequencies][nspec];
    double frequencies[num_frequencies];

    struct Camera *intensityfield;

    init_camera(&intensityfield);

    for (int f = 0; f < num_frequencies; f++) { // For all frequencies...
        for (int s = 0; s < nspec; s++)
            energy_spectrum[f][s] = 0.;
    }

#if (FREQS == FREQLOG)
    for (int f = 0; f < num_frequencies; f++) { // For all frequencies...
        frequencies[f] = FREQ_MIN * pow(10., (double)f / (double)FREQS_PER_DEC);
        fprintf(stderr, "freq = %+.15e\n", frequencies[f]);
    }
#elif (FREQS == FREQFILE)
    FILE *input;
    input = fopen("./frequencies.txt", "r");

    if (input == NULL) {
        fprintf(stderr, "Cannot read frequencies.txt\n");
        exit(1);
    }
    for (int f = 0; f < num_frequencies; f++) {
        fscanf(input, "%lf", &frequencies[f]);
        fprintf(stderr, "freq = %+.15e\n", frequencies[f]);
    }
#endif

    fprintf(stderr, "\nStarting ray tracing\n\n");

#if (SMR)
    prerun_refine(&intensityfield);
#endif

    int block = 0;

    while (block  < tot_blocks) { // block_total
        if (block % (25) == 0)
            fprintf(stderr, "block %d of total %d\n", block, tot_blocks);

        calculate_image_block(&intensityfield[block], frequencies);
#if (AMR)
        if (refine_block(intensityfield[block])) {
            add_block(&intensityfield, block);
        } else {
            block++;
        }
#else
        block++;
#endif
    }

    fprintf(stderr, "\nRay tracing done!\n\n");

    compute_spec(intensityfield, energy_spectrum);

#if (USERSPEC)
    compute_spec_user(intensityfield, energy_spectrum);
#endif
    // WRITE OUTPUT FILES
    /////////////////////

    output_files(intensityfield, energy_spectrum, frequencies);

#if (UNIF)
    write_uniform_camera(intensityfield, frequencies[0], 0);
#endif
    // FREE ALLOCATED POINTERS
    //////////////////////////

    free(intensityfield);

    fprintf(stderr, "\nThat's all folks!\n");

    // END OF PROGRAM
    /////////////////

    return 0;
}
