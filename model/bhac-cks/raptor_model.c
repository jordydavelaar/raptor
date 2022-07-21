
/*
 * raptor_BHAC_AMR_model.c
 *
 */

#include "raptor_model.h"
#include "functions.h"
#include "parameters.h"

void init_model() {
    // Set physical units
    set_units(M_UNIT);

    // Initialize the BHAC AMR GRMHD data
    fprintf(stderr, "\nStarting read in of BHAC GRMHD data...\n");

    init_grmhd_data(GRMHD_FILE);
}

int find_igrid(double x[4], struct block *block_info, double ***Xc) {
    double small = 1e-6;

    for (int igrid = 0; igrid < nleafs; igrid++) {
        if (x[1] + small >= block_info[igrid].lb[0] &&
            x[1] + small <
                block_info[igrid].lb[0] + (block_info[igrid].size[0]) *
                                              block_info[igrid].dxc_block[0] &&
            x[2] + small >= block_info[igrid].lb[1] &&
            x[2] + small <
                block_info[igrid].lb[1] + (block_info[igrid].size[1]) *
                                              block_info[igrid].dxc_block[1] &&
            x[3] + small >= block_info[igrid].lb[2] &&
            x[3] + small <
                block_info[igrid].lb[2] + (block_info[igrid].size[2]) *
                                              block_info[igrid].dxc_block[2]) {
            return igrid;
        }
    }
    return -1;
}

int find_cell(double x[4], struct block *block_info, int igrid, double ***Xc) {

    int i = (int)((x[1] - block_info[igrid].lb[0]) /
                  block_info[igrid].dxc_block[0]);
    int j = (int)((x[2] - block_info[igrid].lb[1]) /
                  block_info[igrid].dxc_block[1]);
    int k = (int)((x[3] - block_info[igrid].lb[2]) /
                  block_info[igrid].dxc_block[2]);

    int cell = i + j * block_info[igrid].size[0] +
               k * block_info[igrid].size[0] * block_info[igrid].size[1];

    return cell;
}

void new_index(int dims, int i, int *new_i, int *new_j, int *new_k, int ip,
               int jp, int kp) {

    if (dims == 1) {
        *new_i = 2 * (ip) + i;
        *new_j = 0;
        *new_k = 0;
    }

    if (dims == 2) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + i / 2;
        *new_k = 0;
    }

    if (dims == 3) {
        *new_i = 2 * (ip) + i % 2;
        *new_j = 2 * (jp) + ((int)(i / 2.)) % 2;
        *new_k = 2 * (kp) + (int)(i / 4.);
    }
}

void read_node(FILE *file_id, int *igrid, int *refine, int ndimini, int level,
               int ind_i, int ind_j, int ind_k) {
    int buffer_i[1], leaf;
    fread(buffer_i, sizeof(int), 1, file_id);
    leaf = buffer_i[0];
    //	fprintf(stderr,"%d\n",leaf);
    if (leaf) {
        (*igrid)++;
        int i = (*igrid); //+(*refine);
        if (i > forest_size) {
            forest = (int *)realloc(forest, i * sizeof(int));
            forest_size++;
        }
        if (i > block_size) {
            block_info =
                (struct block *)realloc(block_info, i * sizeof(struct block));
            block_size++;
        }

        forest[forest_size - 1] = 1;

        block_info[block_size - 1].ind[0] = ind_i;
        block_info[block_size - 1].ind[1] = ind_j;
        block_info[block_size - 1].ind[2] = ind_k;
        block_info[block_size - 1].level = level;
    } else {
        (*refine) += 1;
        int i = (*igrid) + (*refine);
        // printf("%d %d\n",forest_size,i);

        if (i > forest_size) {
            forest = (int *)realloc(forest, i * sizeof(int) * 2);
            forest_size++;
        }
        //  printf("%d %d\n",sizeof(forest),i*sizeof(int));

        forest[forest_size - 1] = 0;
        //  fprintf(stderr,"Blaat\n");
        for (int ich = 0; ich < pow(2, ndimini); ich++) {
            int cind_j, cind_i, cind_k;
            new_index(ndimini, ich, &cind_i, &cind_j, &cind_k, ind_i, ind_j,
                      ind_k);
            read_node(file_id, igrid, refine, ndimini, level + 1, cind_i,
                      cind_j, cind_k);
        }
    }
}

double get_detgamma(double x, double y, double z) {

    double g_dd[4][4];
    double X_u[4];

    X_u[0] = 0;
    X_u[1] = x;
    X_u[2] = y;
    X_u[3] = z;

    metric_dd(X_u, g_dd);

    double detgamma =
        g_dd[1][1] * (g_dd[2][2] * g_dd[3][3] - g_dd[3][2] * g_dd[2][3]) -
        g_dd[1][2] * (g_dd[2][1] * g_dd[3][3] - g_dd[3][1] * g_dd[2][3]) +
        g_dd[1][3] * (g_dd[2][1] * g_dd[3][2] - g_dd[2][2] * g_dd[3][1]);
#if (DEBUG)
    if (isnan(sqrt(detgamma))) {
        double R2 = x * x + y * y + z * z;
        double a2 = a * a;
        double r2 =
            (R2 - a2 + sqrt((R2 - a2) * (R2 - a2) + 4. * a2 * z * z)) * 0.5;
        double r_current = sqrt(r2);

        fprintf(stderr, "isnan detgam %e %e rc %e\n", sqrt(detgamma), detgamma,
                r_current);
        detgamma = 0;
        exit(1);
    }
#endif
    return sqrt(detgamma);
}

void calc_coord_bar(double *X, double *dxc_block, double *Xbar) {
    double coef_1D[3];
    double coef_2D[3][3];
    double coef_3D[3][3][3];

    double is[3], js[3], ks[3];

    coef_1D[0] = 1.;
    coef_1D[1] = 4.;
    coef_1D[2] = 1.;

    is[0] = -1.;
    is[1] = 0.;
    is[2] = 1.;

    js[0] = -1.;
    js[1] = 0.;
    js[2] = 1.;

    ks[0] = -1.;
    ks[1] = 0.;
    ks[2] = 1.;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            coef_2D[i][j] = coef_1D[i] * coef_1D[j];
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                coef_3D[i][j][k] = coef_2D[i][j] * coef_1D[k];
            }
        }
    }

    double norm = 0;
    double xbar = 0;
    double ybar = 0;
    double zbar = 0;
    double detgamma;

    for (int k = 0; k < 3; k++) {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                detgamma = get_detgamma(X[0] + dxc_block[0] / 2. * is[i],
                                        X[1] + dxc_block[1] / 2. * js[j],
                                        X[2] + dxc_block[2] / 2. * ks[k]);
                norm += detgamma * coef_3D[i][j][k];
                xbar += detgamma * coef_3D[i][j][k] *
                        (X[0] + dxc_block[0] / 2. * is[i]);
                ybar += detgamma * coef_3D[i][j][k] *
                        (X[1] + dxc_block[1] / 2. * js[j]);
                zbar += detgamma * coef_3D[i][j][k] *
                        (X[2] + dxc_block[2] / 2. * ks[k]);
                //			fprintf(stderr,"x %e %e
                //%e\n",X[0]+dxc_block[0]/2. * is[i],X[1]+dxc_block[1]/2. *
                // js[j],X[2]+dxc_block[2]/2. * ks[k]);
            }
        }
    }

    norm *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0);
    xbar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;
    ybar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;
    zbar *= (dxc_block[0] / 6.0) * (dxc_block[1] / 6.0) * (dxc_block[2] / 6.0) /
            norm;

    Xbar[0] = xbar;
    Xbar[1] = ybar;
    Xbar[2] = zbar;

//	fprintf(stderr,"Xbar %e %e %e Xgrid %e %e
//%e\n",Xbar[0],Xbar[1],Xbar[2],X[0],X[1],X[2]);
#if (DEBUG)
    if (isnan(Xbar[0])) {
        fprintf(stderr, "isnan in calc_coord_bar %e %e %e %e\n", Xbar[0],
                Xbar[1], Xbar[2], norm);
        Xbar[0] = 0;
        Xbar[1] = 0;
        Xbar[2] = 0;
        exit(1);
    }
#endif
    //	exit(1);
}

void calc_coord(int c, int *nx, int ndimini, double *lb, double *dxc_block,
                double *X) {
    int local_ig[3];

    local_ig[0] = (int)((c % nx[0]));
    local_ig[1] = (int)(fmod((((double)c) / ((double)nx[0])), (double)nx[1]));
    if (ndimini == 3) {
        local_ig[2] = (int)(((double)c) / ((double)nx[0] * nx[1]));
    }

    for (int i = 0; i < ndimini; i++) {
        X[i] = lb[i] + (local_ig[i] + 0.5) * dxc_block[i]; // cell centered
    }

    if (isnan(X[0])) {
        fprintf(stderr, "isnan in calccoord %lf %d %lf %lf\n", lb[0],
                local_ig[0], dxc_block[0], ((double)c) / (nx[0] * nx[1]));
        exit(1);
    }
}

void convert2prim(double prim[8], double **conserved, int c, double X[3],
                  double Xgrid[3], double dxc[3]) {

    double X_u[4];
    double g_dd[4][4], g_uu[4][4];
    X_u[0] = 0;
    X_u[1] = X[0]; // - dxc[0];
    X_u[2] = X[1]; // - dxc[0];
    X_u[3] = X[2]; // - dxc[0];

    double r_current = get_r(X);
    if (r_current < 1.0)
        return;

    metric_dd(X_u, g_dd);
    metric_uu(X_u, g_uu);
    conserved[S2][c] = conserved[S2][c];
    double BS = conserved[S1][c] * conserved[B1][c] +
                conserved[S2][c] * conserved[B2][c] +
                conserved[S3][c] * conserved[B3][c];
    double Bsq = 0;

    double B_d[4];
    double S_u[4];
    S_u[1] = 0;
    S_u[2] = 0;
    S_u[3] = 0;
    B_d[1] = 0;
    B_d[2] = 0;
    B_d[3] = 0;

    double gamma[4][4];
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            gamma[i][j] = g_uu[i][j] + g_uu[0][i] * g_uu[0][j] / (-g_uu[0][0]);
        }
    }

    for (int j = 1; j < 4; j++) {
        for (int i = 1; i < 4; i++) {
            S_u[j] += gamma[i][j] * conserved[S1 + i - 1][c];
            B_d[j] += g_dd[i][j] * conserved[B1 + i - 1][c];
        }
    }

    for (int i = 1; i < 4; i++) {
        //		for(int j=1;j<4;j++){
        Bsq += B_d[i] * conserved[B1 + i - 1][c];
    }
    //	}

#if (DEBUG)
    if (isnan(BS) || isnan(Bsq)) {
        fprintf(stderr, "Bsq %e BS %e\n", Bsq, BS);
        fprintf(stderr, "B %e %e %e\n", conserved[B1][c], conserved[B2][c],
                conserved[B3][c]);
        fprintf(stderr, "V %e %e %e\n", conserved[S1][c], conserved[S2][c],
                conserved[S3][c]);

        LOOP_ij fprintf(stderr, "gij %d %d %e\n", i, j, g_dd[i][j]);
        exit(1);
    }
#endif
    prim[KRHO] = conserved[D][c] / conserved[LFAC][c];
    prim[UU] = (neqpar[0] - 1.) / neqpar[0] *
               (conserved[XI][c] / pow(conserved[LFAC][c], 2.) - prim[KRHO]) /
               (neqpar[0] -
                1.); // need internal energy not pressure so extra 1/(gam-1)
    prim[U1] =
        S_u[1] / (conserved[XI][c] + Bsq) +
        conserved[B1][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));
    prim[U2] =
        S_u[2] / (conserved[XI][c] + Bsq) +
        conserved[B2][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));
    prim[U3] =
        S_u[3] / (conserved[XI][c] + Bsq) +
        conserved[B3][c] * BS / (conserved[XI][c] * (conserved[XI][c] + Bsq));

    prim[B1] = conserved[B1][c];
    prim[B2] = conserved[B2][c];
    prim[B3] = conserved[B3][c];

    if (prim[UU] < 0) {
        fprintf(stderr, "UU %e\n", prim[UU]);
        prim[UU] = (conserved[DS][c] / conserved[D][c]) *
                   pow(prim[KRHO], neqpar[0] - 1) /
                   (neqpar[0] -
                    1.); // need internal energy not pressure so extra 1/(gam-1)
        fprintf(stderr, "UU %e\n", prim[UU]);
        if (prim[UU] < 0) {
            fprintf(stderr, "Entropy UU re_dery failed %e\n", prim[UU]);
            prim[UU] = 1e-14;
            fprintf(stderr, "UU set to low number %e\n", prim[UU]);
        }
    }

    double VdotV = 0;

    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            VdotV += g_dd[i][j] * prim[U1 + i - 1] * prim[U1 + j - 1];
        }
    }
    double gammaf = 1. / sqrt(1 - VdotV);

    if (VdotV > 1.)
        fprintf(stderr, "issues with conserved %e %e %e\n", VdotV, gammaf,
                conserved[LFAC][c]);

#if (DEBUG)

    if (prim[UU] < 0) {
        fprintf(stderr, "U %e gam %e XI %e LFAC %e lor %e RHO %e\n", prim[UU],
                neqpar[0], conserved[XI][c], conserved[LFAC][c], gammaf,
                prim[KRHO]);
        exit(1);
    }

    if (isnan(gammaf)) {
        fprintf(stderr, "VdotV %e lfac %e lor %e\n", VdotV, conserved[LFAC][c],
                gammaf);
        fprintf(stderr, "\n");

        fprintf(stderr, "lor %e vdotv %e lfac %e XI %e Bsq %e BS %e\n", gammaf,
                VdotV, conserved[LFAC][c], conserved[XI][c], Bsq, BS);
        fprintf(stderr, "xi? %e %e\n",
                conserved[LFAC][c] * conserved[LFAC][c] * prim[KRHO] *
                    (1 + neqpar[0] * prim[UU] / prim[KRHO]),
                conserved[XI][c]);
        fprintf(stderr, "rc %e %e\n", r_current, (1. + sqrt(1. - a * a)));
        fprintf(stderr, "Xbar %e %e %e\n", X[0], X[1], X[2]);
        fprintf(stderr, "X %e %e %e\n", Xgrid[0], Xgrid[1], Xgrid[2]);
        fprintf(stderr, "dxc %e %e %e\n", dxc[0], dxc[1], dxc[2]);
        //		exit(1);
    }
    if (isnan(gammaf))
        exit(1);
#endif
}

void init_grmhd_data(char *fname) {

    /*
            nxlone1         = 192
            nxlone2         = 96
            xprobmin1         = 0.19325057145871735
            xprobmax1         = 7.824046010856292
            xprobmin2         = 0.0d0
            xprobmax2         = 0.5d0
    */
    ng[0] = 1;
    ng[1] = 1;
    ng[2] = 1;
    hslope = 1.0;

#if (metric == MKSBHAC || metric == MKSN)
    nxlone[0] = 128;
    nxlone[1] = 48;
    nxlone[2] = 48;

    xprobmax[0] = 8.1117280833;
    xprobmax[1] = M_PI;
    xprobmax[2] = 2. * M_PI;

    xprobmin[0] = 0.17;
    xprobmin[1] = 0.;
    xprobmin[2] = 0.;

#elif (metric == CKS)

    nxlone[0] = 96;
    nxlone[1] = 96;
    nxlone[2] = 384;

    xprobmax[0] = 500;
    xprobmax[1] = 500;
    xprobmax[2] = 2000.;

    xprobmin[0] = -500;
    xprobmin[1] = -500;
    xprobmin[2] = -2000;
#else

    fprintf(stderr, "Not a metric supported by BHAC. Aborting...\n");
    exit(1);

#endif

    startx[1] = xprobmin[0];
    startx[2] = xprobmin[1];
    startx[3] = xprobmin[2];

    stopx[1] = xprobmax[0];
    stopx[2] = xprobmax[1];
    stopx[3] = xprobmax[2];

    double buffer[1];
    unsigned int buffer_i[1];
    FILE *file_id;

    file_id = fopen(fname, "rb"); // r for read, b for binary

    if (file_id < 0) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        fflush(stderr);
        exit(1234);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }

    long int offset;

    int levmaxini, ndirini, nwini, nws, neqparini, it;
    double t;
    fseek(file_id, 0, SEEK_END);
    offset = -36 - 4; // -4; // 7 int, 1 double = 7 * 4 + 1*8 = 56?
    fseek(file_id, offset, SEEK_CUR);

    fread(buffer_i, sizeof(int), 1, file_id);
    nleafs = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    levmaxini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    ndimini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    ndirini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    nwini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    nws = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    neqparini = buffer_i[0];

    fread(buffer_i, sizeof(int), 1, file_id);
    it = buffer_i[0];

    fread(buffer, sizeof(double), 1, file_id);
    t = buffer[0];

    offset = offset - (ndimini * 4 + neqparini * 8);
    fseek(file_id, offset, SEEK_CUR);

    neqpar = (double *)malloc(neqparini * sizeof(double));
    nx = (int *)malloc(ndimini * sizeof(int));

    for (int k = 0; k < ndimini; k++) {
        fread(buffer_i, sizeof(int), 1, file_id);
        nx[k] = buffer_i[0];
    }
    for (int k = 0; k < neqparini; k++) {
        fread(buffer, sizeof(double), 1, file_id);
        neqpar[k] = buffer[0];
    }

    a = neqpar[3];
    Q = 0.0;
    // Q=0.66;//662912607362388;
    int cells = 1;
    for (int k = 0; k < ndimini; k++) {
        cells *= nx[k];
    }

    N1 = nleafs;
    N2 = cells;
    N3 = 1;

    if (nleafs != N1 || cells != N2 || 1 != N3) {
        fprintf(stderr, "wrong N1, N2, N3! %d!=%d %d!=%d %d!=%d", nleafs, N1,
                cells, N2, 1, N3);
    }

    // fprintf(stderr,"cells %d\n",cells);
    long int size_block = cells * (nwini)*8; // size of a block in bytes

    // Read forest
    offset = nleafs * size_block +
             nleafs * (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1) * nws * 8;
    // fprintf(stderr,"%d\n",cells*nwini*nleafs);
    // exit(1);
    fseek(file_id, 0, SEEK_SET);
    fseek(file_id, offset, SEEK_CUR);

    for (int i = 0; i < ndimini; i++) {
        ng[i] = nxlone[i] / nx[i]; // number of blocks in each direction
    }
    if (ndimini < 3)
        ng[2] = 1;
    int igrid = 0;
    int refine = 0;

    block_info = (struct block *)malloc(0);

    forest = (int *)malloc(sizeof(int));
    fprintf(stderr, ".");

    int level = 1;
    for (int k = 0; k < ng[2]; k++) {
        for (int j = 0; j < ng[1]; j++) {
            for (int i = 0; i < ng[0]; i++) {
                //				fprintf(stderr,"%d %d
                //%d\n",i,j,k);
                read_node(file_id, &igrid, &refine, ndimini, level, i, j, k);
            }
        }
    }
    if (nleafs != igrid) {
        fprintf(stderr, "something wrong with grid dimensions\n");
        exit(1);
    }

    double *dx1, *dxc;
    dx1 = (double *)malloc(ndimini * sizeof(double));
    dxc = (double *)malloc(ndimini * sizeof(double));
    // block and cell size on level one.
    for (int i = 0; i < ndimini; i++) {
        dx1[i] = (xprobmax[i] - xprobmin[i]) / ng[i];
        dxc[i] = (xprobmax[i] - xprobmin[i]) / nxlone[i];
    }

    double **values;
    values = (double **)malloc(nwini * sizeof(double *));
    for (int i = 0; i < nwini; i++) {
        values[i] = (double *)malloc(cells * sizeof(double));
    }

    Xgrid = (double ***)malloc(nleafs * sizeof(double **));
    Xbar = (double ***)malloc(nleafs * sizeof(double **));
    for (int j = 0; j < nleafs; j++) {
        Xgrid[j] = (double **)malloc(cells * sizeof(double *));
        Xbar[j] = (double **)malloc(cells * sizeof(double *));
        for (int i = 0; i < cells; i++) {
            Xgrid[j][i] = (double *)malloc(ndimini * sizeof(double));
            Xbar[j][i] = (double *)malloc(ndimini * sizeof(double));
        }
    }

    init_storage();

    fprintf(stderr, ".");

    fseek(file_id, 0, SEEK_SET);

    for (int i = 0; i < nleafs; i++) {
        for (int n = 0; n < ndimini; n++) {
            block_info[i].lb[n] =
                (xprobmin[n] + (block_info[i].ind[n]) * dx1[n] /
                                   pow(2., (double)block_info[i].level - 1.));
#if (DEBUG)
            if (isnan(block_info[i].lb[n])) {
                fprintf(stderr, "NaN %d %lf %d", i,
                        pow(2, block_info[i].level - 1), block_info[i].level);
                exit(1);
            }
#endif
            block_info[i].dxc_block[n] =
                (dxc[n] / (pow(2., (double)block_info[i].level - 1.)));
            block_info[i].size[n] = nx[n];
            //		printf("%e %d %e\n",
            // block_info[i].dxc_block[n],block_info[i].level,dxc[n]);
        }

        for (int nw = 0; nw < nwini; nw++) {
            for (int c = 0; c < cells; c++) {
                fread(buffer, sizeof(double), 1, file_id);
                values[nw][c] = buffer[0];
            }
        }

#pragma omp parallel for shared(values, p) schedule(static, 1)
        for (int c = 0; c < cells; c++) {
            calc_coord(c, nx, ndimini, block_info[i].lb,
                       block_info[i].dxc_block, Xgrid[i][c]);
            double X_u[4] = {0, Xgrid[i][c][0], Xgrid[i][c][1], Xgrid[i][c][2]};
            double r = get_r(X_u);
            if (r > 1.0) {

                calc_coord_bar(Xgrid[i][c], block_info[i].dxc_block,
                               Xbar[i][c]);

#if (DEBUG)
                if (isnan(Xgrid[i][c][0])) {
                    fprintf(stderr, "%d %d", c, i);
                    exit(1);
                }
#endif
                double prim[8];
                convert2prim(prim, values, c, Xbar[i][c], Xgrid[i][c],
                             block_info[i].dxc_block);

                p[KRHO][i][c][0] = prim[KRHO];
                p[UU][i][c][0] = prim[UU];

                p[U1][i][c][0] = prim[U1];
                p[U2][i][c][0] = prim[U2];
                p[U3][i][c][0] = prim[U3];

                p[B1][i][c][0] = prim[B1];
                p[B2][i][c][0] = prim[B2];
                p[B3][i][c][0] = prim[B3];
            }
            if (i == int(nleafs / 2))
                fprintf(stderr, ".");
        }
#pragma omp barrier

        offset = (nx[0] + 1) * (nx[1] + 1) * (nx[2] + 1) * nws * 8;
        fseek(file_id, offset, SEEK_CUR);
    }
    fprintf(stderr, "Done\n");
    free(values);
    free(forest);

    // exit(1);
}

void set_units(double M_unit_) {
    //	double MBH;

    /* set black hole mass */
    /** could be read in from file here,
        along with M_unit and other parameters **/
    //	MBH = 4.e6;

    /** input parameters appropriate to Sgr A* **/
    // double BH_MASS = MBH * MSUN;

    /** from this, calculate units of length, time, mass,
        and derivative units **/
    L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    T_unit = L_unit / SPEED_OF_LIGHT;

    fprintf(stderr, "\nUNITS\n");
    fprintf(stderr, "L,T,M: %g %g %g\n", L_unit, T_unit, M_unit_);

    RHO_unit = M_unit_ / pow(L_unit, 3);
    U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);

    fprintf(stderr, "rho,u,B: %g %g %g\n", RHO_unit, U_unit, B_unit);

    Ne_unit = RHO_unit / (PROTON_MASS + ELECTRON_MASS);
}

void init_storage() {
    int i;
    fprintf(stderr, "\nAllocation memory...\n");
    p = (double ****)malloc(
        NPRIM * sizeof(double ***)); // malloc_rank1(NPRIM, sizeof(double *));
    for (i = 0; i < NPRIM; i++) {
        p[i] = (double ***)malloc(
            (N1 + 1) * sizeof(double **)); // malloc_rank2_cont(N1, N2);
        for (int j = 0; j <= N1; j++) {
            p[i][j] = (double **)malloc((N2 + 1) * sizeof(double *));
            for (int k = 0; k <= N2; k++) {
                p[i][j][k] = (double *)malloc((N3) * sizeof(double));
            }
        }
    }

    values = (double **)malloc(nwini * sizeof(double *));
    for (int j = 0; j < nwini; j++) {
        values[j] = (double *)malloc(cells * sizeof(double));
    }

    return;
}

void coefficients(double X[NDIM], struct block *block_info, int igrid, int c,
                  double del[NDIM]) {

    double block_start[4];
    double block_dx[4];
    int i, j, k;

    block_start[0] =
        block_info[igrid].lb[0] + 0.5 * block_info[igrid].dxc_block[0];
    block_dx[0] = block_info[igrid].dxc_block[0];

    block_start[1] =
        block_info[igrid].lb[1] + 0.5 * block_info[igrid].dxc_block[1];
    block_dx[1] = block_info[igrid].dxc_block[1];

    block_start[2] =
        block_info[igrid].lb[2] + 0.5 * block_info[igrid].dxc_block[2];
    block_dx[2] = block_info[igrid].dxc_block[2];

    i = (int)((X[1] - block_start[0]) / block_dx[0] - 0.5 + 1000) - 1000;
    j = (int)((X[2] - block_start[1]) / block_dx[1] - 0.5 + 1000) - 1000;
    k = (int)((X[3] - block_start[2]) / block_dx[2] + 1000 - 0.5) - 1000;

    if (i < 0) {
        del[1] = 0.;
    } else if (i > nx[0] - 2) {
        //        del[1] = 1.;
        del[1] =
            (X[1] - ((i + 0.5) * block_dx[0] + block_start[0])) / block_dx[0];

    } else {
        del[1] =
            (X[1] - ((i + 0.5) * block_dx[0] + block_start[0])) / block_dx[0];
    }

    if (j < 0) {
        del[2] = 0.;
    } else if (j > nx[1] - 2) {
        //        del[2] = 1.;
        del[2] =
            (X[2] - ((j + 0.5) * block_dx[1] + block_start[1])) / block_dx[1];

    } else {
        del[2] =
            (X[2] - ((j + 0.5) * block_dx[1] + block_start[1])) / block_dx[1];
    }

    if (k < 0) {
        del[3] = 0.;
    } else if (k > nx[2] - 2) {
        //        del[3] = 1.;
        del[3] =
            (X[3] - ((k + 0.5) * block_dx[2] + block_start[2])) / block_dx[2];
    } else {
        del[3] =
            (X[3] - ((k + 0.5) * block_dx[2] + block_start[2])) / block_dx[2];
    }

    return;
}

int compute_c(int i, int j, int k) {
    // computes "c" index for p array given a set of i,j,k indices
    return i + j * nx[0] + k * nx[0] * nx[1];
}

double interp_scalar(double **var, int c, double coeff[4]) {

    double interp;
    int c_ip, c_jp, c_kp;
    int c_i, c_j, c_k = 0;
    double b1, b2, b3, del[NDIM];
    int cindex[2][2][2];

    del[1] = coeff[1];
    del[2] = coeff[2];
    del[3] = coeff[3];
    if (del[1] > 1 || del[2] > 1 || del[3] > 1 || del[1] < 0 || del[2] < 0 ||
        del[3] < 0)
        fprintf(stderr, "del[1] %e \n del[2] %e\n del[3] %e\n", del[1], del[2],
                del[3]);

    c_i = (int)((c % nx[0]));
    c_j = (int)(fmod((((double)c) / ((double)nx[0])), (double)nx[1]));
    if (ndimini == 3) {
        c_k = (int)(((double)c) / ((double)nx[0] * nx[1]));
    }

    c_ip = c_i + 1;
    c_jp = c_j + 1;
    c_kp = c_k + 1;

    if (c_ip >= nx[0]) {
        c_ip = c_i;
        c_i--;
    }

    if (c_jp >= nx[1]) {
        c_jp = c_j;
        c_j--;
    }

    if (c_kp >= nx[2]) {
        c_kp = c_k;
        c_k--;
    }

    b1 = 1. - del[1];
    b2 = 1. - del[2];
    b3 = 1. - del[3];

    cindex[0][0][0] = compute_c(c_i, c_j, c_k);
    cindex[1][0][0] = compute_c(c_ip, c_j, c_k);
    cindex[0][1][0] = compute_c(c_i, c_jp, c_k);
    cindex[0][0][1] = compute_c(c_i, c_j, c_kp);
    cindex[1][1][0] = compute_c(c_ip, c_jp, c_k);
    cindex[1][0][1] = compute_c(c_ip, c_j, c_kp);
    cindex[0][1][1] = compute_c(c_i, c_jp, c_kp);
    cindex[1][1][1] = compute_c(c_ip, c_jp, c_kp);

    /*
        fprintf(stderr,"c %d\n",c);
        fprintf(stderr,"c_i %d c_j %d c_k %d\n",c_i,c_j,c_k);
        fprintf(stderr,"c_ip %d c_jp %d c_kp %d\n",c_ip,c_jp,c_kp);

        fprintf(stderr,"000 %d\n",cindex[0][0][0]);
        fprintf(stderr,"100 %d\n",cindex[1][0][0]);
        fprintf(stderr,"010 %d\n",cindex[0][1][0]);
        fprintf(stderr,"001 %d\n",cindex[0][0][1]);
        fprintf(stderr,"110 %d\n",cindex[1][1][0]);
        fprintf(stderr,"101 %d\n",cindex[1][0][1]);
        fprintf(stderr,"011 %d\n",cindex[0][1][1]);
        fprintf(stderr,"111 %d\n",cindex[1][1][1]);
    */
    interp = var[cindex[0][0][0]][0] * b1 * b2 +
             var[cindex[0][1][0]][0] * b1 * del[2] +
             var[cindex[1][0][0]][0] * del[1] * b2 +
             var[cindex[1][1][0]][0] * del[1] * del[2];

    /* Now interpolate above in x3 */
    interp = b3 * interp + del[3] * (var[cindex[0][0][1]][0] * b1 * b2 +
                                     var[cindex[0][1][1]][0] * b1 * del[2] +
                                     var[cindex[1][0][1]][0] * del[1] * b2 +
                                     var[cindex[1][1][1]][0] * del[1] * del[2]);

    return interp;
}

// Get the fluid parameters in the local co-moving plasma frame.
int get_fluid_params(double X[NDIM], struct GRMHD *modvar) {
    double g_dd[NDIM][NDIM];
    double g_uu[NDIM][NDIM];
    int igrid = (*modvar).igrid_c;
    int i, c;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], V_u[NDIM], VdotV;

#if (metric == MKSBHAC || metric == MKSN)
    X[3] = fmod(X[3], 2 * M_PI);
    X[2] = fmod(X[2], M_PI);
    if (X[3] < 0.)
        X[3] = 2. * M_PI + X[3];
    if (X[2] < 0.)
        X[2] = M_PI + X[2];
#endif

    double r = get_r(X);
    if (r < 1.00)
        return 0;

    double smalll = 1.e-6;
    double small = 0;

    if (X[1] > stopx[1] || X[1] < startx[1] || X[2] < startx[2] ||
        X[2] > stopx[2] || X[3] < startx[3] || X[3] > stopx[3]) {
        return 0;
    }

    if (igrid == -1 || X[1] + small < block_info[igrid].lb[0] ||
        X[1] + small >
            block_info[igrid].lb[0] +
                (block_info[igrid].size[0]) * block_info[igrid].dxc_block[0] ||
        X[2] + small < block_info[igrid].lb[1] ||
        X[2] + small >
            block_info[igrid].lb[1] +
                (block_info[igrid].size[1]) * block_info[igrid].dxc_block[1] ||
        X[3] + small < block_info[igrid].lb[2] ||
        X[3] + small >
            block_info[igrid].lb[2] +
                (block_info[igrid].size[2]) * block_info[igrid].dxc_block[2]) {
        (*modvar).igrid_c = find_igrid(X, block_info, Xgrid);
        igrid = (*modvar).igrid_c;
    }

    if (igrid == -1) {
        fprintf(stderr, "issues with finding igrid, exiting... %e %e %e\n",
                X[1], X[2], X[3]);
        exit(1);
    }

    (*modvar).dx_local = block_info[igrid].dxc_block[0];

    c = find_cell(X, block_info, igrid, Xgrid);

    metric_uu(X, g_uu);

    metric_dd(X, g_dd);

    coefficients(X, block_info, igrid, c, del);

    rho = interp_scalar(p[KRHO][igrid], c, del);
    uu = interp_scalar(p[UU][igrid], c, del);

    (*modvar).n_e = rho * Ne_unit + smalll;

    Bp[1] = interp_scalar(p[B1][igrid], c, del);
    Bp[2] = interp_scalar(p[B2][igrid], c, del);
    Bp[3] = interp_scalar(p[B3][igrid], c, del);
    V_u[1] = interp_scalar(p[U1][igrid], c, del);
    V_u[2] = interp_scalar(p[U2][igrid], c, del);
    V_u[3] = interp_scalar(p[U3][igrid], c, del);

    double gamma_dd[4][4];
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            gamma_dd[i][j] = g_dd[i][j];
        }
    }
    double shift[4];
    for (int j = 1; j < 4; j++) {
        shift[j] = g_uu[0][j] / (-g_uu[0][0]);
    }
    double alpha = 1 / sqrt(-g_uu[0][0]);
    VdotV = 0.;
    (*modvar).U_u[0] = 0.;
    for (int i = 1; i < NDIM; i++) {
        for (int j = 1; j < NDIM; j++) {
            VdotV += gamma_dd[i][j] * V_u[i] * V_u[j];
        }
    }

    if (VdotV > 1.) {
        // fprintf(stderr,
        //         "VdotV too large %e %d %d %e %e %e %e %e %e %e %e %e %e\n",
        //        VdotV, igrid, c, X[1], X[2], X[3], r, V_u[1], V_u[2], V_u[3],
        //       p[U1][igrid][c][0], p[U2][igrid][c][0], p[U3][igrid][c][0]);
        //	issue with normalization, either due to inaccurate
        // interpolation/extrapolation typically only occurs very close to the
        // horizon setting it to a high value lfac=10
        LOOP_i V_u[i] *= sqrt(0.99 / VdotV);
        VdotV = 0.99;
    }

    double lfac = 1 / sqrt(1 - VdotV);

    (*modvar).U_u[0] = lfac / alpha;

    for (int i = 1; i < NDIM; i++) {
        (*modvar).U_u[i] = 0;
        (*modvar).U_u[i] = V_u[i] * lfac - shift[i] * lfac / alpha;
    }

    lower_index(X, (*modvar).U_u, (*modvar).U_d);

    (*modvar).B_u[0] = 0;
    for (i = 1; i < NDIM; i++) {
        for (int l = 1; l < NDIM; l++) {
            (*modvar).B_u[0] +=
                lfac * Bp[i] * (gamma_dd[i][l] * V_u[l]) / alpha;
        }
    }

    for (int i = 1; i < NDIM; i++) {
        (*modvar).B_u[i] = 0;
        (*modvar).B_u[i] =
            (Bp[i] + alpha * (*modvar).B_u[0] * (*modvar).U_u[i]) / lfac;
    }

    lower_index(X, (*modvar).B_u, (*modvar).B_d);

    // magnetic field
    double Bsq = fabs((*modvar).B_u[0] * (*modvar).B_d[0] +
                      (*modvar).B_u[1] * (*modvar).B_d[1] +
                      (*modvar).B_u[2] * (*modvar).B_d[2] +
                      (*modvar).B_u[3] * (*modvar).B_d[3]) +
                 smalll;
    (*modvar).B = sqrt(Bsq) * B_unit;

#if (DEBUG)
    if (isnan(Bsq) || isnan((*modvar).B)) {
        double R2 = Xgrid[igrid][c][1] * Xgrid[igrid][c][1] +
                    Xgrid[igrid][c][2] * Xgrid[igrid][c][2] +
                    Xgrid[igrid][c][3] * Xgrid[igrid][c][3];
        double a2 = a * a;
        double r2 = (R2 - a2 +
                     sqrt((R2 - a2) * (R2 - a2) +
                          4. * a2 * Xgrid[igrid][c][3] * Xgrid[igrid][c][3])) *
                    0.5;
        fprintf(stderr, "B isnan r %e rmin %e\n", sqrt(r2), CUTOFF_INNER);
        fprintf(stderr, "B isnan X %e %e %e %e\n", X[0], X[1], X[2], X[3]);
        fprintf(stderr, "B isnan Xgrid %e %e %e\n", Xgrid[igrid][c][0],
                Xgrid[igrid][c][1], Xgrid[igrid][c][2]);
        fprintf(stderr, "B isnan Bsq %e B_u %e %e %e %e U_u %e %e %e %e\n", Bsq,
                (*modvar).B_u[0], (*modvar).B_u[1], (*modvar).B_u[2],
                (*modvar).B_u[3], (*modvar).U_u[0], (*modvar).U_u[1],
                (*modvar).U_u[2], (*modvar).U_u[3]);
        fprintf(stderr, "B isnan Bp %e %e %e\n", Bp[1], Bp[2], Bp[3]);
        fprintf(stderr, "B isnan V_u %e %e %e\n", V_u[1], V_u[2], V_u[3]);
        fprintf(stderr, "B isnan VdotV %e\n", VdotV);
        fprintf(stderr, "B isnan lapse %e\n", sqrt(-g_uu[0][0]));
        fprintf(stderr, "B isnan shift %e %e %e\n", g_uu[0][1], g_uu[0][2],
                g_uu[0][3]);
        exit(1);
    }
#endif
    double gam = neqpar[0];

    double beta_trans = 1.0;

    (*modvar).beta = uu * (gam - 1.) / (0.5 * (Bsq + smalll));

    double b2 = pow(((*modvar).beta / beta_trans), 2.);

    (*modvar).sigma = Bsq / (rho + smalll);

    (*modvar).sigma_min = 1.0;
    double Rhigh = R_HIGH;
    double Rlow = R_LOW;

    double trat = Rhigh * b2 / (1. + b2) + Rlow / (1. + b2);

    Thetae_unit = (gam - 1.) * (MPoME) / (trat + 1);

    (*modvar).theta_e = (uu / rho) * Thetae_unit;

    if ((Bsq / (rho + 1e-20) > 5.) || r > 2000. ||
        (*modvar).theta_e > 100.) { // excludes all spine emmission
        (*modvar).n_e = 0;
        return 0;
    }

    return 1;
}
