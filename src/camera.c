/*
 * Radboud Polarized Integrator
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka
 *
 */
#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////

int num_blocks, tot_blocks;
double BLOCK_SIZE_X, BLOCK_SIZE_Y;
int max_level;

// FUNCTIONS
////////////

// Initializes the camera
void init_camera(struct Camera **intensityfield) {
    int x, y;

    num_blocks = IMG_HEIGHT / num_pixels_1d;
    int num_blocks2 = IMG_WIDTH / num_pixels_1d;
    if (num_blocks2 != num_blocks) {
        fprintf(stderr, "should be squared. EXITING\n");
        exit(1);
    }
    tot_blocks = num_blocks * num_blocks2;
    (*intensityfield) = malloc((tot_blocks) * sizeof(struct Camera));
    for (int block = 0; block < tot_blocks; block++) {
        (*intensityfield)[block].level = 1;

        x = (int)block / (double)num_blocks;
        y = block % num_blocks;
        (*intensityfield)[block].ind[0] = x;
        (*intensityfield)[block].ind[1] = y;

        get_impact_params(intensityfield, block);
    }
}

// Iniatlized the impact parameters for every pixel in a block
void get_impact_params(struct Camera **intensityfield, int block) {
    int xpixel, ypixel;
    double stepx, stepy, d_x, d_y;

    BLOCK_SIZE_X = CAM_SIZE_X / (pow(2, (*intensityfield)[block].level - 1) *
                                 (double)(num_blocks));
    BLOCK_SIZE_Y = CAM_SIZE_Y / (pow(2, (*intensityfield)[block].level - 1) *
                                 (double)(num_blocks));

    d_x = BLOCK_SIZE_X * R_GRAV / SOURCE_DIST; // angular size of block in cm
    d_y = BLOCK_SIZE_Y * R_GRAV / SOURCE_DIST; // angular size of block in cm

    (*intensityfield)[block].lcorner[0] =
        -CAM_SIZE_X * 0.5 + (*intensityfield)[block].ind[0] * BLOCK_SIZE_X;
    (*intensityfield)[block].lcorner[1] =
        -CAM_SIZE_Y * 0.5 + (*intensityfield)[block].ind[1] * BLOCK_SIZE_Y;

    (*intensityfield)[block].dx[0] = (d_x / (double)num_pixels_1d);
    (*intensityfield)[block].dx[1] = (d_y / (double)num_pixels_1d);

    for (int pixel = 0; pixel < tot_pixels; pixel++) {
        xpixel = (int)pixel / (double)num_pixels_1d;
        ypixel = pixel % num_pixels_1d;
        stepx = BLOCK_SIZE_X / ((double)num_pixels_1d);
        stepy = BLOCK_SIZE_Y / ((double)num_pixels_1d);
        (*intensityfield)[block].alpha[pixel] =
            (xpixel + 0.5) * stepx + (*intensityfield)[block].lcorner[0];
        (*intensityfield)[block].beta[pixel] =
            (ypixel + 0.5) * stepy + (*intensityfield)[block].lcorner[1];

        for (int f = 0; f < num_frequencies; f++) {
            for (int s = 0; s < 4; s++) {
                (*intensityfield)[block].IQUV[pixel][f][s] = 0;
            }
            (*intensityfield)[block].tau[pixel][f] = 0;
            (*intensityfield)[block].tauF[pixel][f] = 0;
        }
    }
}

// For new Adapative Grid Block it compute the new indices.
void new_cindex(int child, int *new_i, int *new_j, int ip, int jp) {
    *new_i = 2 * (ip) + child % 2;
    *new_j = 2 * (jp) + child / 2;
}

// Shift array so that there is space for the new block
void shift_camera_array(struct Camera **intensityfield, int current_block) {
    for (int block = tot_blocks - 1; block > current_block + 3; block--) {
        (*intensityfield)[block] = (*intensityfield)[block - 3];
    }
}

// Splits the original block in a new set of four blocks in the camera struct
void add_block(struct Camera **intensityfield, int current_block) {
    int cind_i, cind_j;

    int ind_i = (*intensityfield)[current_block].ind[0];
    int ind_j = (*intensityfield)[current_block].ind[1];
    tot_blocks += 3;

    (*intensityfield) =
        realloc((*intensityfield), (tot_blocks) * sizeof(struct Camera));

    shift_camera_array(intensityfield, current_block);
    // compute new indices
    int new_level = (*intensityfield)[current_block].level + 1;
    for (int i = 0; i < 4; i++) {
        new_cindex(i, &cind_i, &cind_j, ind_i, ind_j);
        (*intensityfield)[current_block + i].ind[0] = cind_i;
        (*intensityfield)[current_block + i].ind[1] = cind_j;
        (*intensityfield)[current_block + i].level = new_level;

        get_impact_params(intensityfield, current_block + i);
    }
}

// Checks if a refinement criterion is met, returns 1 if that is the case
int refine_block(struct Camera intensity) {
    int pixel1, pixel2, pixel3;
    double gradI_x, gradI_y;
    double gradImax = -1e100;

    for (int xpixel = 0; xpixel < num_pixels_1d - 1; xpixel++) {
        for (int ypixel = 0; ypixel < num_pixels_1d - 1; ypixel++) {
            for (int freq = 0; freq < num_frequencies; freq++) {
                pixel1 = ypixel + xpixel * num_pixels_1d;
                pixel2 = ypixel + 1 + xpixel * num_pixels_1d;
                pixel3 = ypixel + (xpixel + 1) * num_pixels_1d;

                gradI_y = fabs(intensity.IQUV[pixel2][freq][0] -
                               intensity.IQUV[pixel1][freq][0]) /
                          (intensity.IQUV[pixel1][freq][0] + 1e-40);
                gradI_x = fabs(intensity.IQUV[pixel3][freq][0] -
                               intensity.IQUV[pixel1][freq][0]) /
                          (intensity.IQUV[pixel1][freq][0] + 1e-40);

                if (gradI_x > gradImax)
                    gradImax = gradI_x;
                if (gradI_y > gradImax)
                    gradImax = gradI_y;
            }
        }
    }

    if (gradImax > 0.25 && intensity.level < max_level)
        return 1;
    else
        return 0;
}

// Static Camera Grid, checks if a block should be refined before ray tracing
// begins
int refine_init_block(struct Camera intensity) {

    double radius_5 = 20;
    double radius_4 = 30;
    double radius_3 = 40;
    double radius_2 = 60;

    BLOCK_SIZE_X =
        CAM_SIZE_X / (pow(2, intensity.level - 1) * (double)(num_blocks));
    BLOCK_SIZE_Y =
        CAM_SIZE_Y / (pow(2, intensity.level - 1) * (double)(num_blocks));

    double lcorner_x = intensity.lcorner[0];
    double lcorner_y = intensity.lcorner[1];
    double ucorner_x = intensity.lcorner[0] + BLOCK_SIZE_X;
    double ucorner_y = intensity.lcorner[1] + BLOCK_SIZE_Y;

    double rl = sqrt(lcorner_x * lcorner_x + lcorner_y * lcorner_y);
    double ru = sqrt(ucorner_x * ucorner_x + ucorner_y * ucorner_y);

    double rmin = fmin(rl, ru);

    bool bool_5 = radius_5 > rmin;
    bool bool_4 = radius_4 > rmin && radius_5 < rmin;
    bool bool_3 = radius_3 > rmin && radius_4 < rmin;
    bool bool_2 = radius_2 > rmin && radius_3 < rmin;

    if (bool_5 && intensity.level < 5 && max_level > 4)
        return 1;
    if (bool_4 && intensity.level < 4 && max_level > 3)
        return 1;
    if (bool_3 && intensity.level < 3 && max_level > 2)
        return 1;
    if (bool_2 && intensity.level < 2 && max_level > 1)
        return 1;
    else
        return 0;
}

// Goes over all blocks before ray tracing and adds new block if refinement
// criterion is met
void prerun_refine(struct Camera **intensityfield) {
    int block = 0;
    while (block < tot_blocks) {
        if (refine_init_block((*intensityfield)[block])) {
            add_block((intensityfield), block);
        } else {
            block++;
        }
    }
}

// Initialzies a single pixel, assigns wave vector to it.
void init_pixel(double alpha, double beta, double t, double photon_u[8]) {
#if (LINEAR_IMPACT_CAM)
    double stepx = CAM_SIZE_X / (double)IMG_WIDTH;
    double stepy = CAM_SIZE_Y / (double)IMG_HEIGHT;
    alpha = -CAM_SIZE_X * 0.5 + (x + 0.5) * stepx;
    beta = -CAM_SIZE_Y * 0.5 + (y + 0.5) * stepy;

    initialize_photon(alpha, beta, photon_u, t_init);
#endif

#if (LOG_TETRADS_CAM || LINEAR_TETRADS_CAM)
    initialize_photon_tetrads(alpha, beta, photon_u, t, Xcam, Ucam);
#endif
}

// Given an impact parameter, finds the corresponding block that it is in
int find_block(double x[2], struct Camera *intensityfield) {
    double small = 1e-6;
    double dx[2];

    for (int block = 0; block < tot_blocks; block++) {
        dx[0] = intensityfield[block].dx[0] * SOURCE_DIST / R_GRAV;
        dx[1] = intensityfield[block].dx[1] * SOURCE_DIST / R_GRAV;
        // fprintf(stderr,"x  = %lf, %lf \n", x[0],x[1]);
        // fprintf(stderr,"dx = %lf, %lf \n", dx[0], dx[1]);
        // fprintf(stderr,"nulpunt2 = %lf \n",
        // intensityfield[block].lcorner[0]);
        if (x[0] + small >= intensityfield[block].lcorner[0] &&
            x[0] + small <
                num_pixels_1d * dx[0] + intensityfield[block].lcorner[0] &&
            x[1] + small >= intensityfield[block].lcorner[1] &&
            x[1] + small <
                num_pixels_1d * dx[1] + intensityfield[block].lcorner[1]) {
            return block;
        }
    }

    return -1;
}

// Given a impact parameter, finds the pixels it is closest to
int find_pixel(double x[2], struct Camera *intensityfield, int block) {
    double dx[2];
    dx[0] = intensityfield[block].dx[0] * SOURCE_DIST / R_GRAV;
    dx[1] = intensityfield[block].dx[1] * SOURCE_DIST / R_GRAV;

    int i = (int)((x[0] - intensityfield[block].lcorner[0]) / dx[0] - 0.5);
    int j = (int)((x[1] - intensityfield[block].lcorner[1]) / dx[1] - 0.5);

    int pixel = j + i * num_pixels_1d;

    return pixel;
}
