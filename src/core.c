/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////

char GRMHD_FILE[256];

double MBH, M_UNIT, TIME_INIT, INCLINATION;
double R_HIGH, R_LOW;
double FREQS_PER_DEC, FREQ_MIN, FREQ_MAX;

double SOURCE_DIST; // Distance to M87 (cm); for Sgr A* use (2.47e22)

int IMG_WIDTH, IMG_HEIGHT;
double CAM_SIZE_X, CAM_SIZE_Y;
double STEPSIZE;

// FUNCTIONS
////////////

// Read model parameters from model.in
void read_model(char *argv[]) {
    char temp[100], temp2[100];
    FILE *input;
    char inputfile[100];

    sscanf(argv[1], "%s", inputfile);
    fprintf(stdout, "\nUsing model parameter file %s\n", inputfile);

    input = fopen(inputfile, "r");
    if (input == NULL) {
        fprintf(stderr, "Can't read file %s! Aborting", inputfile);
        exit(1);
    }

    // Model parameters
    fscanf(input, "%s %s %lf", temp, temp2, &MBH);
    fscanf(input, "%s %s %lf", temp, temp2, &SOURCE_DIST);
    fscanf(input, "%s %s %lf", temp, temp2, &M_UNIT);
    fscanf(input, "%s %s %lf", temp, temp2, &R_LOW);
    fscanf(input, "%s %s %lf", temp, temp2, &R_HIGH);
    fscanf(input, "%s %s %lf", temp, temp2, &INCLINATION);

    // Observer parameters
    fscanf(input, "%s %s %d", temp, temp2, &IMG_WIDTH);
    fscanf(input, "%s %s %d", temp, temp2, &IMG_HEIGHT);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_X);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_Y);

    fscanf(input, "%s %s %lf", temp, temp2, &FREQS_PER_DEC);
    fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MIN);
    fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE);
    fscanf(input, "%s %s %d", temp, temp2, &max_level);

    // Second argument: GRMHD file
    sscanf(argv[2], "%s", GRMHD_FILE);
    sscanf(argv[3], "%lf", &TIME_INIT);

    fprintf(stderr, "\nModel parameters:\n\n");
    fprintf(stderr, "MBH \t\t= %g Msun\n", MBH);
    fprintf(stderr, "DISTANCE \t= %g kpc\n", SOURCE_DIST);
    fprintf(stderr, "M_UNIT \t\t= %g grams\n", M_UNIT);
    fprintf(stderr, "R_LOW \t\t= %g \n", R_LOW);
    fprintf(stderr, "R_HIGH \t\t= %g \n", R_HIGH);
    fprintf(stderr, "INCLINATION \t= %g deg\n", INCLINATION);

    fprintf(stderr, "METRIC \t\t= ");
#if (metric == MKSBHAC)
    fprintf(stderr, "MKS BHAC\n");
#elif (metric == CKS)
    fprintf(stderr, "CKS BHAC\n");
#elif (metric == MKSHARM)
    fprintf(stderr, "MKS HARM3D\n");
#endif

    fprintf(stderr, "\nObserver parameters:\n\n");
    fprintf(stderr, "IMG_WIDTH \t= %d \n", IMG_WIDTH);
    fprintf(stderr, "IMG_HEIGHT \t= %d \n", IMG_HEIGHT);
    fprintf(stderr, "CAM_SIZE_X \t= %g GM/c2\n", CAM_SIZE_X);
    fprintf(stderr, "CAM_SIZE_Y \t= %g GM/c2\n", CAM_SIZE_Y);
    fprintf(stderr, "FREQS_PER_DEC \t= %lf \n", FREQS_PER_DEC);
    fprintf(stderr, "FREQ_MIN \t= %g Hz\n", FREQ_MIN);
    fprintf(stderr, "STEPSIZE \t= %g \n", STEPSIZE);

    // to cgs units
    MBH *= MSUN;
    SOURCE_DIST *= KPCTOCM;

    fclose(input);
}

// For a single block this function will iterate over the pixels and call
// geodesic integrations as well as radiation transfer
void calculate_image_block(struct Camera *intensityfield,
                           double frequencies[num_frequencies], int block) {

#pragma omp parallel for shared(frequencies, intensityfield, p)                \
    schedule(static, 1)
    for (int pixel = 0; pixel < tot_pixels; pixel++) {
        int steps = 0;

        double *lightpath2 = malloc(9 * max_steps * sizeof(double));

#if (POL)
        double f_x = 0.;
        double f_y = 0.;
        double p = 0.;
#endif
        // INTEGRATE THIS PIXEL'S GEODESIC
        integrate_geodesic((*intensityfield).alpha[pixel],
                           (*intensityfield).beta[pixel], lightpath2, &steps,
                           CUTOFF_INNER);
        // PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
#if (POL)
        for (int f = 0; f < num_frequencies; f++) {

            radiative_transfer_polarized(lightpath2, steps, frequencies[f],
                                         (*intensityfield).IQUV[pixel][f],
                                         &(*intensityfield).tau[pixel][f],
                                         &(*intensityfield).tauF[pixel][f],
                                         &(*intensityfield).pdf[pixel][f],
                                         (*intensityfield).avg[pixel][f],
                                         block, pixel);
        }

#else
        radiative_transfer_unpolarized(lightpath2, steps, frequencies,
                                       (*intensityfield).IQUV[pixel],
                                       &(*intensityfield).tau[pixel]);
        for (int f = 0; f < num_frequencies; f++) {
            (*intensityfield).IQUV[pixel][f][0] *= pow(frequencies[f], 3.);
        }
#endif
        free(lightpath2);
    }
#pragma omp barrier
}

// Functions that computes a spectrum at every frequency
// by integrating over the image struct
void compute_spec(struct Camera *intensityfield,
                  double energy_spectrum[num_frequencies][nspec]) {
    double dA, S_I, S_Q, S_U, S_V;

    for (int block = 0; block < tot_blocks; block++) {
        dA = (intensityfield)[block].dx[0] * (intensityfield)[block].dx[1];
        for (int pixel = 0; pixel < tot_pixels; pixel++) {
            for (int freq = 0; freq < num_frequencies; freq++) {
#if (POL)

                S_I = (intensityfield)[block].IQUV[pixel][freq][0];
                S_Q = (intensityfield)[block].IQUV[pixel][freq][1];
                S_U = (intensityfield)[block].IQUV[pixel][freq][2];
                S_V = (intensityfield)[block].IQUV[pixel][freq][3];

                // Stokes I
                energy_spectrum[freq][0] += S_I * dA;

                // Stokes Q
                energy_spectrum[freq][1] += S_Q * dA;

                // Stokes U
                energy_spectrum[freq][2] += S_U * dA;

                // stokes V
                energy_spectrum[freq][3] += S_V * dA;

#else
                energy_spectrum[freq][0] +=
                    (intensityfield)[block].IQUV[pixel][freq][0] * dA;
#endif
            }
        }
    }
}
