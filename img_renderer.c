/*
 * Radboud Polarized Integrator v1.0
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Moscibrodzka, Ziri Younsi, Jordy Davelaar
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "functions.h"
#include "constants.h"
#include "parameters.h"

int main(int argc, char *argv[]){

    // Initialize RNG
    init_RCARRY(13242834);


    // INPUT FILE
    /////////////

    char inputfile[100];
// model to read
    sscanf(argv[1], "%s", inputfile);

    FILE *input;
    input = fopen(inputfile, "r");
    if (input == NULL){
        printf ("Cannot read input file");
        return 1;
    }

    char temp[100], temp2[100];

//    read_in_table("symphony_pure_thermal.txt");

    // Model parameters
    fscanf(input, "%s %s %lf", temp, temp2, &MBH);
    fscanf(input, "%s %s %lf", temp, temp2, &M_UNIT);
    fscanf(input, "%s %s %d",  temp, temp2, &ABSORPTION);
    fscanf(input, "%s %s %s",  temp, temp2, TEMP_MODEL);
    fscanf(input, "%s %s %d",  temp, temp2, &SPHERICAL_ACC);

    // Observer parameters
    fscanf(input, "%s %s %d",  temp, temp2, &IMG_WIDTH);
    fscanf(input, "%s %s %d",  temp, temp2, &IMG_HEIGHT);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_X);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_Y);
    fscanf(input, "%s %s %d",  temp, temp2, &FREQS_PER_DEC);
    fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MIN);
    fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MAX);
    fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE);

    // Second argument: GRMHD file
    sscanf(argv[2], "%s", GRMHD_FILE);

    // 3th and 4th arguments: inclination, t_0
    sscanf(argv[3], "%lf", &INCLINATION);
    sscanf(argv[4], "%lf", &TIME_INIT);

    printf("Model parameters:\n");
    printf("MBH \t\t= %g \n", MBH);
    printf("M_UNIT \t\t= %g \n", M_UNIT);
    printf("ABSORPTION \t= %d \n", ABSORPTION);
    printf("TEMP_MODEL \t= %s \n", TEMP_MODEL);
    printf("SPHERICAL_ACC \t= %d \n\n", SPHERICAL_ACC);

    printf("Observer parameters:\n");
    printf("IMG_WIDTH \t= %d \n", IMG_WIDTH);
    printf("IMG_HEIGHT \t= %d \n", IMG_HEIGHT);
    printf("CAM_SIZE_X \t= %g \n", CAM_SIZE_X);
    printf("CAM_SIZE_Y \t= %g \n", CAM_SIZE_Y);
    printf("FREQS_PER_DEC \t= %d \n", FREQS_PER_DEC);
    printf("FREQ_MIN \t= %g \n", FREQ_MIN);
    printf("FREQ_MAX \t= %g \n", FREQ_MAX);
    printf("INCLINATION \t= %g \n", INCLINATION);
    printf("STEPSIZE \t= %g \n", STEPSIZE);
    fclose (input);

    // INITIALIZE MODEL
    ///////////////////

    // Initialize HARM2D grmhd model
    // Note: this sets the black hole spin 'a'
    init_model();

    // Set constants such as R_ISCO, JANSKY_FACTOR
    // These depend on the black hole spin
    set_constants();

printf("\n L_UNIT = %+.15e", L_unit);
printf("\n T_UNIT = %+.15e", T_unit);

    // INITIALIZE VARIABLES
    ///////////////////////

    // Set the "camera time" (initial coordinate time of rays)
//    double t_init = 0.; // MUST BE UPDATED
//    sscanf(argv[1], "%lf", &t_init);

    // Field of intensity values for output
    double lambdafield[IMG_WIDTH * IMG_HEIGHT];
    int    steps;

    // Stepsize for constructing the impact parameters alpha, beta
    double stepx = CAM_SIZE_X / (double) IMG_WIDTH;
    double stepy = CAM_SIZE_Y / (double) IMG_HEIGHT;
    double photon_u[8], alpha, beta;
    int    x, y, f;

    // MAIN PROGRAM LOOP
    ////////////////////

// Change this to:
// initialize_ray(x, y, photon_u);
// img[y * IMG_WIDTH + x] = integrate_backward(photon_u);

    int num_indices = FREQS_PER_DEC * (int) (log10(FREQ_MAX) - log10(FREQ_MIN)) + 1;
    printf("\nNumber of frequencies to compute: %d\n", num_indices);
    double energy_spectrum[num_indices];
    double frequencies[num_indices];
    double intensityfield[num_indices][IMG_WIDTH * IMG_HEIGHT];
    double f_x_field[IMG_WIDTH * IMG_HEIGHT];
    double f_y_field[IMG_WIDTH * IMG_HEIGHT];
    double p_field[IMG_WIDTH * IMG_HEIGHT];
    double I_field[IMG_WIDTH * IMG_HEIGHT];
    double Q_field[IMG_WIDTH * IMG_HEIGHT];
    double U_field[IMG_WIDTH * IMG_HEIGHT];
    double V_field[IMG_WIDTH * IMG_HEIGHT];

    for(f = 0; f < num_indices; f++){ // For all frequencies...
        frequencies[f] = FREQ_MIN * pow(10., (double) f / (double) FREQS_PER_DEC);
        energy_spectrum[f] = 0.;
        printf("freq = %+.15e\n", frequencies[f]);
    }

    for(x = 0; x < IMG_WIDTH; x++){ // For all pixel columns...
        #pragma omp parallel for default(none) private(f,steps,alpha,beta,photon_u) shared(num_indices,energy_spectrum,frequencies,intensityfield,f_x_field,f_y_field,p_field,lambdafield,I_field, Q_field, U_field, V_field,x,stepx,stepy,CUTOFF_INNER, IMG_WIDTH, IMG_HEIGHT, CAM_SIZE_X, CAM_SIZE_Y) schedule(static,1)
        for(y = 0; y < IMG_HEIGHT; y++){ // For all pixel rows (distributed over threads)...
		fprintf(stderr,"%d %d\n",x,y);
            double *lightpath2 = malloc(9 * max_steps * sizeof(double));
            double *IQUV = malloc(4 * sizeof(double));
            double f_x, f_y, p;

            // Compute impact parameters for this pixel
            alpha = -CAM_SIZE_X * 0.5 + (x + 0.5) * stepx;
            beta  = -CAM_SIZE_Y * 0.5 + (y + 0.5) * stepy;

            // INTEGRATE THIS PIXEL'S GEODESIC

            intensityfield[0][y * IMG_WIDTH + x] = integrate_geodesic(alpha, beta, photon_u, lightpath2, &steps, CUTOFF_INNER, &f_x, &f_y, &p, IQUV);
            f_x_field[y * IMG_WIDTH + x] = f_x;
            f_y_field[y * IMG_WIDTH + x] = f_y;
            p_field[y * IMG_WIDTH + x] = p;

            // PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
            for(f = 0; f < num_indices; f++){
//                intensityfield[f][y * IMG_WIDTH + x] = radiative_transfer(lightpath2, steps, frequencies[f]);
                energy_spectrum[f] += intensityfield[f][y * IMG_WIDTH + x];
            }

            I_field[y * IMG_WIDTH + x] = IQUV[0];
            Q_field[y * IMG_WIDTH + x] = IQUV[1];
            U_field[y * IMG_WIDTH + x] = IQUV[2];
            V_field[y * IMG_WIDTH + x] = IQUV[3];

            free(lightpath2);
	    free(IQUV);
        }
        #pragma omp barrier
    }

    // WRITE OUTPUT FILES
    /////////////////////

    // We open ONE spectrum file and multiple image files (one per frequency)
//    FILE *spectrum    = fopen("output/spectrum.dat", "w");

    for(f = 0; f < num_indices;f++){ // For all frequencies...
        // Create filenames, open files
        char dat_filename[256] = "";
        char vtk_filename[256] = "";
        sprintf(dat_filename,"output/img_data_%e_IQUV.dat",frequencies[f]);
        sprintf(vtk_filename,"output/img_data_%e.vtk",frequencies[f]);
        FILE *imgfile     = fopen(dat_filename, "w");
  //      FILE *fp          = fopen(vtk_filename, "w");

        // Write image data to file
//        write_image(imgfile, intensityfield[f], f_x_field, f_y_field, p_field, JANSKY_FACTOR);
        write_image_IQUV(imgfile, I_field, Q_field, U_field, V_field, JANSKY_FACTOR);
//        write_VTK_image(fp, intensityfield[f], lambdafield, JANSKY_FACTOR);

        // Close image files
        fclose(imgfile);
  //      fclose(fp);

        // Update spectrum file
//        fprintf(spectrum, "%+.15e\t%+.15e\n", frequencies[f], JANSKY_FACTOR * energy_spectrum[f]);
    }

//    fclose(spectrum);

    // END OF PROGRAM
    /////////////////

    return 0;
}
