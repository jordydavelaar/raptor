/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "constants.h"
#include "parameters.h"
#include <math.h>
#include <stdio.h>

void write_image(FILE *imgfile, double *intensityfield, double scalefactor) {
    int i, j;

    // Write image to output file
    for (i = 0; i < IMG_WIDTH; i++) {
        for (j = 0; j < IMG_HEIGHT; j++) {
            fprintf(imgfile, "%d\t%d\t%+.15e\n", i, j,
                    scalefactor * intensityfield[i + j * IMG_WIDTH]);
        }
    }
}

void write_image_polarized(FILE *imgfile, double *intensityfield,
                           double *f_x_field, double *f_y_field,
                           double *p_field, double scalefactor) {
    int i, j;

    // Write image to output file
    for (i = 0; i < IMG_WIDTH; i++) {
        for (j = 0; j < IMG_HEIGHT; j++) {
            fprintf(imgfile, "%d\t%d\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", i, j,
                    scalefactor * intensityfield[i + j * IMG_WIDTH],
                    f_x_field[i + j * IMG_WIDTH], f_y_field[i + j * IMG_WIDTH],
                    p_field[i + j * IMG_WIDTH]);
        }
    }
}

void write_image_IQUV(FILE *imgfile, double *Ifield, double *Qfield,
                      double *Ufield, double *Vfield, double scalefactor) {
    int i, j;

    // Write image to output file
    for (i = 0; i < IMG_WIDTH; i++) {
        for (j = 0; j < IMG_HEIGHT; j++) {
            fprintf(imgfile, "%d\t%d\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n", i, j,
                    scalefactor * Ifield[i + j * IMG_WIDTH],
                    scalefactor * Qfield[i + j * IMG_WIDTH],
                    scalefactor * Ufield[i + j * IMG_WIDTH],
                    scalefactor * Vfield[i + j * IMG_WIDTH]);
        }
    }
}

void write_VTK_image(FILE *fp, double *intensityfield, double *lambdafield,
                     double scalefactor) {
    int i, j;
    double stepx = CAM_SIZE_X / (double)IMG_WIDTH;
    double stepy = CAM_SIZE_Y / (double)IMG_HEIGHT;
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Image Simulation Result\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", IMG_WIDTH, IMG_HEIGHT, 1);
    fprintf(fp, "POINTS %d float\n", IMG_WIDTH * IMG_HEIGHT);
    for (j = 0; j < IMG_WIDTH; j++)
        for (i = 0; i < IMG_HEIGHT; i++) {
            double xx = -CAM_SIZE_X * 0.5 + (i + 0.5) * stepx;
            double yy = -CAM_SIZE_Y * 0.5 + (j + 0.5) * stepy;
            fprintf(fp, "%g %g %g\n", xx, yy, 0.0);
        }
    fprintf(fp, "\nPOINT_DATA %d\n", IMG_WIDTH * IMG_HEIGHT);
    fprintf(fp, "SCALARS Intensity float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    double flux = 0.0;
    for (j = 0; j < IMG_WIDTH; j++)
        for (i = 0; i < IMG_HEIGHT; i++) {
            flux += scalefactor * intensityfield[i + j * IMG_WIDTH];
            fprintf(fp, "%+.15e\n",
                    scalefactor * intensityfield[i + j * IMG_WIDTH]);
        }
    fprintf(fp, "SCALARS lambda float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (j = 0; j < IMG_WIDTH; j++)
        for (i = 0; i < IMG_WIDTH; i++) {
            fprintf(fp, "%+.15e\n", lambdafield[i + j * IMG_WIDTH]);
        }
    fprintf(stdout, "Integrated flux density = %.5e\n", flux);
}
