/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "constants.h"
#include "parameters.h"
#include <math.h>
#include <stdio.h>

void read_model(char *argv[]) {
    char inputfile[100];
    // model to read
    sscanf(argv[1], "%s", inputfile);

    FILE *input;
    input = fopen(inputfile, "r");
    if (input == NULL) {
        printf("Cannot read input file");
        //return 1;
    }

    char temp[100], temp2[100];

    //    read_in_table("symphony_pure_thermal.txt");

    // Model parameters
    fscanf(input, "%s %s %lf", temp, temp2, &MBH);
    fscanf(input, "%s %s %lf", temp, temp2, &M_UNIT);
    fscanf(input, "%s %s %d", temp, temp2, &ABSORPTION);
    fscanf(input, "%s %s %s", temp, temp2, TEMP_MODEL);
    fscanf(input, "%s %s %d", temp, temp2, &SPHERICAL_ACC);

    // Observer parameters
    fscanf(input, "%s %s %d", temp, temp2, &IMG_WIDTH);
    fscanf(input, "%s %s %d", temp, temp2, &IMG_HEIGHT);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_X);
    fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_Y);
    fscanf(input, "%s %s %d", temp, temp2, &FREQS_PER_DEC);
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
    fclose(input);
}
