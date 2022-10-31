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

// FUNCTIONS
////////////

void output_files(struct Camera *intensityfield,
                  double energy_spectrum[num_frequencies][nspec],
                  double frequencies[num_frequencies]) {
    struct stat st = {0};
    char spec_folder[64] = "output";

    if (stat(spec_folder, &st) == -1) {
        mkdir(spec_folder, 0700);
    }

#if (SPECFILE)
    char spec_filename[256] = "";
    sprintf(spec_filename, "%s/spectrum_%d_%.02lf.dat", spec_folder,
            (int)TIME_INIT, INCLINATION);
    FILE *specfile = fopen(spec_filename, "w");
#endif

    for (int f = 0; f < num_frequencies; f++) { // For all frequencies...
                                                // Create filenames, open files

#if (IMGFILE)
        char hdf5_filename[512] = "";
        sprintf(hdf5_filename, "%s/img_data_%d.h5", spec_folder,
                (int)TIME_INIT);
        write_image_hdf5(hdf5_filename, intensityfield, frequencies,
                         JANSKY_FACTOR);
#endif

        fprintf(stderr, "Frequency %.5e Hz Integrated flux density = %.5e Jy\n",
                frequencies[f], JANSKY_FACTOR * energy_spectrum[f][0]);
        fprintf(stderr, "Frequency %.5e Hz Bol Luminosity = %.5e ergs/s\n",
                frequencies[f],
                frequencies[f] * energy_spectrum[f][0] * JANSKY_FACTOR * 4 *
                    M_PI * SOURCE_DIST * SOURCE_DIST / 1e23);

#if (SPECFILE)
#if (POL)
        fprintf(specfile, "%+.15e", frequencies[f]);
        for (int s = 0; s < nspec; s++) {
            fprintf(specfile, "\t%+.15e",
                    JANSKY_FACTOR * energy_spectrum[f][s]);
        }
        fprintf(specfile, "\n");
#else
        fprintf(specfile, "%+.15e\t%+.15e\n", frequencies[f],
                JANSKY_FACTOR * energy_spectrum[f][0]);
#endif
#endif
    }
#if (SPECFILE)
    fclose(specfile);
#endif
}

void write_image_hdf5(char *hdf5_filename, struct Camera *data,
                      double *frequencies, double factor) {

    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    herr_t status;
    double dA;

    file_id = H5Fcreate(hdf5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    double buffer[tot_blocks][tot_pixels];

    dims[1] = tot_pixels;
    dims[0] = tot_blocks;

    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];

        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1;
                buffer[block][pixel] =
                    data[block].IQUV[pixel][freq][0] * factor * dA;
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "I%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }

    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];
        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1;
                buffer[block][pixel] = data[block].tau[pixel][freq];
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "tau%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }

#if (POL)
    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];

        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1;
                buffer[block][pixel] =
                    data[block].IQUV[pixel][freq][1] * factor * dA;
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "Q%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }

    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];

        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1; // data[block].dx[0] * data[block].dx[1];
                buffer[block][pixel] =
                    data[block].IQUV[pixel][freq][2] * factor * dA;
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "U%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }

    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];

        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1;
                buffer[block][pixel] =
                    data[block].IQUV[pixel][freq][3] * factor * dA;
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "V%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }

    for (int freq = 0; freq < num_frequencies; freq++) {
        char dataset[200];
        for (int block = 0; block < tot_blocks; block++) {
            for (int pixel = 0; pixel < tot_pixels; pixel++) {
                dA = 1;
                buffer[block][pixel] = data[block].tauF[pixel][freq];
            }
        }

        dataspace_id = H5Screate_simple(2, dims, NULL);

        sprintf(dataset, "tauF%e", frequencies[freq]);
        dataset_id =
            H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 buffer);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
    }
#endif

    char dataset[200];

    dataspace_id = H5Screate_simple(2, dims, NULL);

    for (int block = 0; block < tot_blocks; block++) {
        for (int pixel = 0; pixel < tot_pixels; pixel++) {
            buffer[block][pixel] = data[block].alpha[pixel];
        }
    }
    sprintf(dataset, "alpha");
    dataset_id = H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             buffer);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(2, dims, NULL);

    for (int block = 0; block < tot_blocks; block++) {
        for (int pixel = 0; pixel < tot_pixels; pixel++) {
            buffer[block][pixel] = data[block].beta[pixel];
        }
    }
    sprintf(dataset, "beta");
    dataset_id = H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             buffer);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);

    status = H5Fclose(file_id);
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

// Outputs the ACG camera struct with a uniform resolution
void write_uniform_camera(struct Camera *intensityfield, double frequency,
                          int freq) {
    int uniform_size = IMG_WIDTH * pow(2, max_level - 1);
    double uniform_dx = CAM_SIZE_X / (double)uniform_size;
    double x[2];
    char spec_folder[64] = "output";
    char uniform_filename[256] = "";
    sprintf(uniform_filename, "%s/uniform_img_%.02e_%d.dat", spec_folder,
            frequency, (int)TIME_INIT);

    FILE *uniformfile = fopen(uniform_filename, "w");

    double UNIT_FACTOR = JANSKY_FACTOR * uniform_dx * uniform_dx * R_GRAV *
                         R_GRAV / SOURCE_DIST / SOURCE_DIST;
    double arcsec_factor = 206265.0 * R_GRAV / SOURCE_DIST;

    for (int i = 0; i < uniform_size; i++) {
        for (int j = 0; j < uniform_size; j++) {
            x[0] = -CAM_SIZE_X / 2. + (i + 0.5) * uniform_dx;
            x[1] = -CAM_SIZE_X / 2. + (j + 0.5) * uniform_dx;
            int block = find_block(x, intensityfield);
            if (block == -1) {
                fprintf(stderr, "camera block not found!\n");
                exit(1);
            }
            int pixel = find_pixel(x, intensityfield, block);
            fprintf(uniformfile,
                    "%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\t%+.15e\n",
                    x[0] * arcsec_factor, x[1] * arcsec_factor,
                    intensityfield[block].IQUV[pixel][freq][0] * UNIT_FACTOR,
                    intensityfield[block].IQUV[pixel][freq][1] * UNIT_FACTOR,
                    intensityfield[block].IQUV[pixel][freq][2] * UNIT_FACTOR,
                    intensityfield[block].IQUV[pixel][freq][3] * UNIT_FACTOR,
                    intensityfield[block].tau[pixel][freq],
                    intensityfield[block].tauF[pixel][freq]);
        }
    }
    fclose(uniformfile);
}
