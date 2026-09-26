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

// SMR on the critical curve (off unless definitions.h enables it)
#ifndef SMR_CRITCURVE
#define SMR_CRITCURVE 0
#endif
#ifndef SMR_MAX_LEVEL
#define SMR_MAX_LEVEL max_level // may exceed max_level, the AMR limit
#endif
#ifndef SMR_BUFFER
#define SMR_BUFFER 1.0 // refinement band half-width, in block widths
#endif

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
        (*intensityfield)[block].norder[pixel] = 0.;
        (*intensityfield)[block].norder_mino[pixel] = 0.;
        (*intensityfield)[block].ncross[pixel] = 0;
        (*intensityfield)[block].kappa1[pixel] = 0.;
        (*intensityfield)[block].kappa2[pixel] = 0.;
        (*intensityfield)[block].dchi_grav[pixel] = 0.;
        (*intensityfield)[block].geo_dlam[pixel] = 0.;
        (*intensityfield)[block].geo_deta[pixel] = 0.;
        (*intensityfield)[block].geo_thpole[pixel] = 0.;

        for (int f = 0; f < num_frequencies; f++) {
            for (int s = 0; s < 4; s++) {
                (*intensityfield)[block].IQUV[pixel][f][s] = 0;
            }
            (*intensityfield)[block].tau[pixel][f] = 0;
            (*intensityfield)[block].tauF[pixel][f] = 0;
            (*intensityfield)[block].fnorm_dev[pixel][f] = 0;
            (*intensityfield)[block].fnorm_th[pixel][f] = 0;
            (*intensityfield)[block].fnorm_cam[pixel][f] = 0;
            (*intensityfield)[block].wp_dchi[pixel][f] = 0;
            (*intensityfield)[block].wp_dpsi[pixel][f] = 0;
        }
    }
}

// For new Adapative Grid Block it compute the new indices.
void new_cindex(int child, int *new_i, int *new_j, int ip, int jp) {
    *new_i = 2 * (ip) + child % 2;
    *new_j = 2 * (jp) + child / 2;
}

// Splits the original block in a new set of four blocks in the camera struct.
// The first child takes the parent's slot (so the caller re-processes it at
// the same index) and the other three are appended at the end of the array,
// where the main loop reaches them later. This is O(1) per split; inserting
// all four in place meant shifting every later block (O(N^2) overall, hours
// of memmove at ~1e5 blocks). Block order is irrelevant for the output:
// every block carries its own impact parameters.
void add_block(struct Camera **intensityfield, int current_block) {
    int cind_i, cind_j;

    int ind_i = (*intensityfield)[current_block].ind[0];
    int ind_j = (*intensityfield)[current_block].ind[1];
    int new_level = (*intensityfield)[current_block].level + 1;
    int first_new = tot_blocks;
    tot_blocks += 3;

    (*intensityfield) =
        realloc((*intensityfield), (tot_blocks) * sizeof(struct Camera));
    if ((*intensityfield) == NULL) {
        fprintf(stderr, "add_block: cannot allocate %d blocks\n", tot_blocks);
        exit(1);
    }

    // compute new indices
    for (int i = 0; i < 4; i++) {
        int b = (i == 0) ? current_block : first_new + i - 1;
        new_cindex(i, &cind_i, &cind_j, ind_i, ind_j);
        (*intensityfield)[b].ind[0] = cind_i;
        (*intensityfield)[b].ind[1] = cind_j;
        (*intensityfield)[b].level = new_level;

        get_impact_params(intensityfield, b);
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

#if (SMR_CRITCURVE)
// Critical curve (image of the unstable spherical photon orbits for a distant
// observer; Bardeen 1973, Gralla & Lupsasca 2020) as a closed polyline in
// camera coordinates (alpha, beta), with the conventions of
// initialize_photon() in metric.c (lambda = -alpha sin i, p_theta = beta).
#define CC_NHALF 4096
static double cc_x[2 * CC_NHALF + 1], cc_y[2 * CC_NHALF + 1];
static double cc_box[4]; // xmin, xmax, ymin, ymax

// lambda(r) and beta^2(r) of the spherical photon orbit at radius r, a > 0
static double cc_lam_b2(double r, double aa, double th, double *b2) {
    double Del = r * r - 2. * r + aa * aa;
    double lam = aa + r / aa * (r - 2. * Del / (r - 1.));
    double eta = r * r * r / (aa * aa) * (4. * Del / ((r - 1.) * (r - 1.)) - r);
    double cot = cos(th) / sin(th);
    *b2 = eta + aa * aa * cos(th) * cos(th) - lam * lam * cot * cot;
    return lam;
}

static void critical_curve_init(void) {
    double th = INCLINATION / 180. * M_PI;
    double s = a < 0. ? -1. : 1.;
    double aa = fabs(a);
    int n = 2 * CC_NHALF + 1;

    if (aa < 1e-6 || sin(th) < 1e-6) {
        // Circle: Schwarzschild (radius sqrt 27), or face-on observer, where
        // only the orbit with lambda = 0 is seen, at radius sqrt(eta + a^2)
        double rad = sqrt(27.);
        if (aa >= 1e-6) {
            double lo = 1., hi = 4., b2; // lambda(r) decreases through 0
            for (int k = 0; k < 100; k++) {
                double mid = 0.5 * (lo + hi);
                double Del = mid * mid - 2. * mid + aa * aa;
                if (aa + mid / aa * (mid - 2. * Del / (mid - 1.)) < 0.)
                    hi = mid;
                else
                    lo = mid;
            }
            double r = 0.5 * (lo + hi), Del = r * r - 2. * r + aa * aa;
            b2 = r * r * r / (aa * aa) * (4. * Del / ((r - 1.) * (r - 1.)) - r);
            rad = sqrt(b2 + aa * aa);
        }
        for (int k = 0; k < n; k++) {
            double phi = 2. * M_PI * k / (n - 1);
            cc_x[k] = rad * cos(phi);
            cc_y[k] = rad * sin(phi);
        }
    } else {
        // beta^2 >= 0 on [r1, r2] inside the photon shell [r_pro, r_retro];
        // single-humped there, so bisect for both roots from the peak
        double rpro = 2. * (1. + cos(2. / 3. * acos(-aa)));
        double rret = 2. * (1. + cos(2. / 3. * acos(aa)));
        double b2, b2max = -1e100, rpk = rpro;
        for (int k = 0; k <= 2000; k++) {
            double r = rpro + (rret - rpro) * k / 2000.;
            cc_lam_b2(r, aa, th, &b2);
            if (b2 > b2max) {
                b2max = b2;
                rpk = r;
            }
        }
        double ends[2], outer[2] = {rpro, rret};
        for (int e = 0; e < 2; e++) {
            double lo = outer[e], hi = rpk; // b2(lo) < 0 <= b2(hi)
            for (int k = 0; k < 100; k++) {
                double mid = 0.5 * (lo + hi);
                cc_lam_b2(mid, aa, th, &b2);
                if (b2 < 0.)
                    lo = mid;
                else
                    hi = mid;
            }
            ends[e] = hi;
        }
        // cosine spacing clusters points at the roots, where beta ~ sqrt
        for (int k = 0; k < CC_NHALF; k++) {
            double r = ends[0] + (ends[1] - ends[0]) * 0.5 *
                                     (1. - cos(M_PI * k / (CC_NHALF - 1.)));
            double lam = cc_lam_b2(r, aa, th, &b2);
            double beta = sqrt(fmax(b2, 0.));
            cc_x[k] = cc_x[2 * CC_NHALF - 1 - k] = -s * lam / sin(th);
            cc_y[k] = beta;
            cc_y[2 * CC_NHALF - 1 - k] = -beta;
        }
        cc_x[n - 1] = cc_x[0];
        cc_y[n - 1] = cc_y[0];
    }

    cc_box[0] = cc_box[2] = 1e100;
    cc_box[1] = cc_box[3] = -1e100;
    for (int k = 0; k < n; k++) {
        cc_box[0] = fmin(cc_box[0], cc_x[k]);
        cc_box[1] = fmax(cc_box[1], cc_x[k]);
        cc_box[2] = fmin(cc_box[2], cc_y[k]);
        cc_box[3] = fmax(cc_box[3], cc_y[k]);
    }
    fprintf(stderr,
            "SMR: critical curve for a = %g, i = %g deg: alpha in [%g, %g], "
            "beta in [%g, %g]\n",
            a, INCLINATION, cc_box[0], cc_box[1], cc_box[2], cc_box[3]);
}

// 1 if the point (x, y) lies within distance d of the critical curve
static int near_critical_curve(double x, double y, double d) {
    double bx = fmax(fmax(cc_box[0] - x, x - cc_box[1]), 0.);
    double by = fmax(fmax(cc_box[2] - y, y - cc_box[3]), 0.);
    if (bx * bx + by * by >= d * d)
        return 0;
    for (int k = 0; k < 2 * CC_NHALF; k++) {
        double sx = cc_x[k + 1] - cc_x[k], sy = cc_y[k + 1] - cc_y[k];
        double px = x - cc_x[k], py = y - cc_y[k];
        double l2 = sx * sx + sy * sy;
        double t = l2 > 0. ? fmin(fmax((px * sx + py * sy) / l2, 0.), 1.) : 0.;
        double dx = px - t * sx, dy = py - t * sy;
        if (dx * dx + dy * dy < d * d)
            return 1;
    }
    return 0;
}
#endif

// Static Camera Grid, checks if a block should be refined before ray tracing
// begins
int refine_init_block(struct Camera intensity) {

    double block_size_x =
        CAM_SIZE_X / (pow(2, intensity.level - 1) * (double)(num_blocks));
    double block_size_y =
        CAM_SIZE_Y / (pow(2, intensity.level - 1) * (double)(num_blocks));

#if (SMR_CRITCURVE)
    // Refine every block whose centre lies within SMR_BUFFER block widths
    // (plus the half diagonal) of the critical curve. Subrings approach the
    // curve geometrically, so this makes the pixel size proportional to the
    // distance from the curve down to level SMR_MAX_LEVEL.
    double b = fmax(block_size_x, block_size_y);
    double xc = intensity.lcorner[0] + 0.5 * block_size_x;
    double yc = intensity.lcorner[1] + 0.5 * block_size_y;
    return intensity.level < SMR_MAX_LEVEL &&
           near_critical_curve(xc, yc, (SMR_BUFFER + 0.5 * sqrt(2.)) * b);
#else
    double radius_5 = 20;
    double radius_4 = 30;
    double radius_3 = 40;
    double radius_2 = 60;

    double lcorner_x = intensity.lcorner[0];
    double lcorner_y = intensity.lcorner[1];
    double ucorner_x = intensity.lcorner[0] + block_size_x;
    double ucorner_y = intensity.lcorner[1] + block_size_y;

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
#endif
}

// Goes over all blocks before ray tracing and splits every block that meets
// the refinement criterion, one level per pass. Children replace their parent
// in place (same block order as add_block), but each pass copies the array
// once instead of shifting it per split, so this is O(N) per level rather
// than O(N^2).
void prerun_refine(struct Camera **intensityfield) {
#if (SMR_CRITCURVE)
    critical_curve_init();
#endif
    for (;;) {
        int *flag = malloc(tot_blocks * sizeof(int));
        int nref = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : nref)
        for (int block = 0; block < tot_blocks; block++) {
            flag[block] = refine_init_block((*intensityfield)[block]);
            nref += flag[block];
        }
        if (nref == 0) {
            free(flag);
            break;
        }

        int new_tot = tot_blocks + 3 * nref;
        struct Camera *refined = malloc(new_tot * sizeof(struct Camera));
        if (refined == NULL) {
            fprintf(stderr, "prerun_refine: cannot allocate %d blocks\n",
                    new_tot);
            exit(1);
        }
        int j = 0, cind_i, cind_j;
        for (int block = 0; block < tot_blocks; block++) {
            if (!flag[block]) {
                refined[j++] = (*intensityfield)[block];
                continue;
            }
            int ind_i = (*intensityfield)[block].ind[0];
            int ind_j = (*intensityfield)[block].ind[1];
            int new_level = (*intensityfield)[block].level + 1;
            for (int i = 0; i < 4; i++, j++) {
                new_cindex(i, &cind_i, &cind_j, ind_i, ind_j);
                refined[j].ind[0] = cind_i;
                refined[j].ind[1] = cind_j;
                refined[j].level = new_level;
                get_impact_params(&refined, j);
            }
        }
        free(flag);
        free(*intensityfield);
        *intensityfield = refined;
        tot_blocks = new_tot;
        fprintf(stderr, "SMR: split %d blocks, now %d blocks\n", nref,
                tot_blocks);
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
