/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 * A list of all RAPTOR functions.
 */

#include "definitions.h"
#include "model_definitions.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// CORE.C
/////////

void read_model(char *argv[]);

// GRMATH.C
///////////

double get_r(double X_u[4]);

// Lowers the index of the contravariant vector V_u, storing the results in
// a covariant one (V_d), based on the metric at position X_u
void lower_index(double X_u[4], double V_u[4], double V_d[4]);

// Lowers two indices of a rank (2, 0) tensor
void lower_two_indices(double N_uu[4][4], double N_dd[4][4], double X_u[4]);

// Lowers the index of a contravariant vector V_u in BL coordinates.
void BL_lower_index(double X_u[4], double V_u[4], double V_d[4]);

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
void raise_index(double X_u[4], double V_d[4], double V_u[4]);

// Raises the index of the covariant vector V_d, storing the results in a
// contravariant one (V_u), based on the metric at position X_u
// Needed for CKS coordinates
void raise_index_KS(double X_u[4], double V_d[4], double V_u[4]);

// Adjusts y[4] = U_u[0] so that y describes a lightray/null geodesic
void normalize_null(double X_u[4], double U_u[4]);

// Returns the norm of U_u, which is the scalar g_dd[a][b] * U_u[a] * U_u[b]
double four_velocity_norm(double X_u[4], double U_u[4]);

// Returns the inner product of vectors A and B, i.e. A_u B_d
double inner_product(double *X_u, double *A_u, double *B_u);

// Transform a contravariant vector from BL to KS coordinates
void BL_to_KS_u(double *BLphoton_u, double *KSphoton_u);

// Transform a contravariant vector from KS to BL coordinates
void KS_to_BL_u(double *KSphoton_u, double *BLphoton_u);

// Convert KS to CKS coordinates
void KS_to_CKS(double *X_KS_u, double *X_CKS_u);

// Convert CKS to KS coordinates
void CKS_to_KS(double *X_CKS_u, double *X_KS_u);

// Transform a contravariant vector from KS to CKS coordinates
void KS_to_CKS_u(double *KScoords, double *CKScoords);

// Return the photon frequency in the co-moving frame of the plasma
double freq_in_plasma_frame(double Uplasma_u[4], double k_d[4]);

// Angle between k_u and B_u in the plasma frame
double pitch_angle(double *X_u, double *k_u, double *B_u, double *Uplasma_u);

// void f_tetrad_to_stokes(double Iinv, double Iinv_pol, double complex
// f_tetrad_u[], double complex S_A[4]);

// void stokes_to_f_tetrad(double complex S_A[], double *Iinv, double *Iinv_pol,
// double complex f_tetrad_u[4]);

// void construct_U_vector( double X_u[], double U_u[]);

// INTEGRATOR.C
///////////////

// Returns an appropriate stepsize dlambda, which depends on position & velocity
double stepsize(double X_u[4], double U_u[4]);

// Updates the vector y (containing position/velocity) by one RK4 step.
void rk4_step(double *y, void (*f)(double *, double *), double dt);

// Updates the vector y (containing position/velocity) by one RK2 step.
void rk2_step(double *y, void (*f)(double *, double *), double dt);

// Updates the vector y (containing position/velocity) by one Verlet step.
void verlet_step(double *y, void (*f)(double *, double *), double dl);

// The function to be used by the integrator - it solves the geodesic eqn.
// y contains the 4-position and the 4-velocity for one lightray/particle
void f_geodesic(double *y, double *fvector);

// Integrate the null geodesic specified by alpha and beta, store results
// in lightpath
void integrate_geodesic(double alpha, double beta, double *lightpath,
                        int *steps, double cutoff_inner);

void radiative_transfer_polarized(double *lightpath, int steps,
                                  double frequency, double *f_x, double *f_y,
                                  double *p, int PRINT_POLAR, double *IQUV,
                                  double *tau, double *tauF);

double radiative_transfer_unpolarized(double *lightpath, int steps,
                                      double *frequency,
                                      double IQUV[num_frequencies][4],
                                      double tau[num_frequencies]);
// METRIC.C
///////////

// Computes the metric at location X
void metric_dd(double X_u[4], double g_dd[4][4]);

// Computes the inverse metric at location X
void metric_uu(double X_u[4], double g_uu[4][4]);

// Computes the inverse metric at location X
void metric_KS_uu(double X_u[4], double g_uu[4][4]);

// Computes the Christoffel symbols at location X numerically (general metric)
void connection_num_udd(double X_u[4], double gamma_udd[4][4][4]);

// Computes the Christoffel symbols at location X based on an exact metric
void connection_udd(double X_u[4], double gamma_udd[4][4][4]);

// This function initializes a single 'superphoton' or light ray.
void initialize_photon(double alpha, double beta, double k_u[4], double t_init);

// Transformation functions
double Xg2_approx_rand(double Xr2);

double Ug2_approx_rand(double Ur2, double Xg2);

// TETRAD.C
///////////

double determ(double matrix[][4], int n);

void create_tetrad(double X_u[], double k_u[], double U_u[],
                   double tetrad_u[][4]);

void create_observer_tetrad(double X_u[], double k_u[], double U_u[],
                            double b_u[], double tetrad_u[][4]);

double tetrad_identity_eta(double X_u[4], double tetrad_u[4][4], int a, int b);

double tetrad_identity_g(double tetrad_u[][4], int mu, int nu);

double tetrad_identity_sum_latin(double tetrad_u[4][4], double tetrad_d[4][4],
                                 int mu, int nu);

double tetrad_identity_sum_greek(double tetrad_u[4][4], double tetrad_d[4][4],
                                 int a, int b);

void create_tetrad_d(double X_u[], double tetrad_u[][4], double tetrad_d[][4]);

double check_tetrad_compact(double X_u[], double tetrad_u[][4]);

void check_tetrad_identities(double X_u[], double tetrad_u[][4]);

// EMISSION.C
/////////////

void evaluate_coeffs_user(double *jI, double *jQ, double *jU, double *jV,
                          double *rQ, double *rU, double *rV, double *aI,
                          double *aQ, double *aU, double *aV, double nu_p,
                          struct GRMHD modvar, double pitch_ang);

void evaluate_coeffs_single(double *jI, double *jQ, double *jU, double *jV,
                            double *rQ, double *rU, double *rV, double *aI,
                            double *aQ, double *aU, double *aV, double nu_p,
                            struct GRMHD modvar, double pitch_ang);

// Return emission coefficient j_nu for kappa distribution function
double emission_coeff_kappa_FIT(double nu, double Ne, double Thetae, double B,
                                double theta);

// Return absorption coefficient for kappa distribution function
double absorption_coeff_kappa_FIT(double nu, double Ne, double Thetae, double B,
                                  double theta);

// Return emission coefficient j_nu for thermal synchrotron radiation
double emission_coeff_THSYNCH(double B_, double theta, double THETA_e_,
                              double nu_plasma, double n_e);

// Return emission coefficient for angle-averaged thermal synchrotron radiation
double emission_coeff_THSYNCHAV(double B_, double THETA_e_, double nu_plasma,
                                double n_e);

// Return emission coefficient for thermal free-free radiation
double emission_coeff_FFTHERMAL(double nu, double n_e, double T);

// Simple emissivity model: orbiting Gaussian hotspot (see Dexter 2009)
double emissivity_hotspot(double *X_u);

// Simple emissivity model: thin disk line emission (see Dexter 2009)
double emissivity_thindisk(double *X_u);

// Return absorption coefficient a_nu
double absorption_coeff_TH(double j_nu, double nu, double THETA_e);

// Planck function
double planck_function(double nu, double THETA_e);

// Perform radiative transfer along the ray stored in "lightpath"
double radiative_transfer(double *lightpath, int steps, double frequency);

// Backward transfer
double backward_transfer(double alpha, double beta, double *photon_u,
                         int *steps);

// POL_EMISSION.C
/////////////////

double j_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B);
double j_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B);
double j_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B);

double a_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_I_thermal);
double a_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_Q_thermal);
double a_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_V_thermal);

double j_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);
double j_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);
double j_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);

double a_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);
double a_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);
double a_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B);

double rho_V_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B);
double rho_V_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B);
double rho_Q_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B);
double rho_Q_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B);

double j_I(double theta_e, double n_e, double nu, double B, double theta_B);
double j_Q(double theta_e, double n_e, double nu, double B, double theta_B);
double j_V(double theta_e, double n_e, double nu, double B, double theta_B);

double a_I(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_I_thermal);
double a_Q(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_Q_thermal);
double a_V(double theta_e, double n_e, double nu, double B, double theta_B,
           double j_V_thermal);

double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B);
double rho_V(double theta_e, double n_e, double nu, double B, double theta_B);

// UTILITIES.C
//////////////

void set_constants();

// IO.C
///////

// Write the array "intensityfield" (scaled by "scalefactor") to the file
// "imgfile"
void write_image(FILE *imgfile, double *intensityfield, double scalefactor);

void write_image_polarized(FILE *imgfile, double *intensityfield,
                           double *f_x_field, double *f_y_field,
                           double *p_field, double scalefactor);

void write_image_IQUV(FILE *imgfile, double *Ifield, double *Qfield,
                      double *Ufield, double *Vfield, double scalefactor);

// Write the arrays "intensityfield" (scaled by "scalefactor") and "lambdafield"
// to a VTK file
void write_VTK_image(FILE *imgfile, double *intensityfield, double *lambdafield,
                     double scalefactor);

// RAPTOR_HARM_MODEL.C
//////////////////////

// See grmonty paper by Dolence et al.
// HARM model internal utilities
void init_model();

void set_units(double);

void init_grmhd_data(char *fname);

void init_storage();

void Xtoijk(double *X, int *i, int *j, int *k, double *del);

// void get_fluid_params(double X[4], double *Ne, double *Thetae, double *B,
//                      double *B_u, double Ucon[4], int *IN_VOLUME);
int get_fluid_params(double X[NDIM], struct GRMHD *modvar);
// IO

void compute_spec(struct Camera *intensity,
                  double energy_spectrum[num_frequencies][nspec]);

void compute_spec_user(struct Camera *intensity,
                       double energy_spectrum[num_frequencies][nspec]);

// Create output files (image, spectrum, etc.)
void output_files(struct Camera *intesityfield,
                  double spectrum[num_frequencies][nspec],
                  double frequencies[num_frequencies]);

void write_image_hdf5(char *hdf5_filename, struct Camera *data,
                      double *frequencies, double factor);

void write_uniform_camera(struct Camera *intensityfield, double frequency,
                          int freq);
// Integrate null geodesics, perform radiative transfer calculations, and
// compute the image.
void calculate_image_block(struct Camera *intensityfield,
                           double frequencies[num_frequencies]);
/// CAMERA.C
void init_camera(struct Camera **intensityfield);

void add_block(struct Camera **intensityfield, int current_block);

int refine_block();

void prerun_refine(struct Camera **intensityfield);

void get_impact_params(struct Camera **intensityfield, int block);

int find_block(double x[2], struct Camera *intensityfield);

int find_pixel(double x[2], struct Camera *intensityfield, int block);

#endif // FUNCTIONS_H
