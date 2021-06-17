/*
 * Radboud Polarized Integrator official V1.0
 * Copyright 2014-2020 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Monika Mościbrodzka, Ziri Younsi
 *
 * This program integrates the equation of motion for a particle in a
 * gravitational field, defined by the metric included by the user.
 *
 * Indices 0, 1, 2, 3 correspond to t, r, theta, phi (Schwarzschild/Kerr).
 *
 * (Null) geodesics are parametrized by an (affine) parameter called lambda.
 *
 * Sign convention: (-,+,+,+)
 *
 * Indices are labeled: "u" (up) means contravariant index
 *                      "d" (down) means covariant index
 *
 * Examples: U_u[alpha], U_d[alpha], gamma_udd[mu][alpha][beta]
 *
 * Definition of the "y_u vector" which describes a ray (position/velocity):
 * y[0] = X_u[0] // Position
 * y[1] = X_u[1]
 * y[2] = X_u[2]
 * y[3] = X_u[3]
 * y[4] = U_u[0] // Wavevector
 * y[5] = U_u[1]
 * y[6] = U_u[2]
 * y[7] = U_u[3]
 */

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "constants.h"
#include "functions.h"
#include "parameters.h"

int main(int argc, char *argv[]) {
  // Initialize variables
  ///////////////////////

  double E_init, E;
  double L_init, L;
  double Q_init, Q;
  double norm, r_current, t_current;
  double t_final = 1.e5;
  int i;

  double lambda;           // The affine parameter, or proper time for particles
  double dlambda_adaptive; // Affine parameter variable

  // Set the "camera time" (initial coordinate time of rays)
  double t_init;
  sscanf(argv[1], "%lf", &t_init);

  // Compute horizon radius for integration cutoff
  double Rh = (1. + sqrt(1. - a * a));
  double cutoff_inner =
      Rh * (1. + horizon_marg); // Cutoff outside or inside BH EH

  // Field of intensity values for output
  // static arrays because we do not want to call malloc functions within the
  // loop stack variables will work faster two global, shared arrays
  double k_u[4]; // Covariant and contravariant wavevector
  double X_u[4]; // current position
  int steps;

  // Stepsize for constructing the impact parameters alpha, beta
  double photon_u[8], alpha, beta;

  // Output file containing the geodesic
  FILE *GRfile = fopen("GRoutput.dat", "w");
  fprintf(GRfile, "#Lambda \t\tt \t\t\tr \t\t\ttheta \t\t\tphi \t\t\tUt \
                     Ur \t\t\tUtheta \t\t\tUphi \t\t\tE \t\t\tL \t\t\tQ \
                     \tNorm \n");

  double gammatest[4][4][4];
  double gammatestnum[4][4][4];
  double Xtest_u[4] = {100., 1.8134, 0.642, 1.9877};
  connection_udd(Xtest_u, gammatest);
  connection_num_udd(Xtest_u, gammatestnum);

  // Initialize the photon for this pixel
  alpha = 5.3; //-0.001;
  beta = 0.0;  // 4.752;
  initialize_photon(alpha, beta, photon_u, t_init);
  // initialize_photon_perspective(alpha, beta, photon_u, t_init);
  photon_u[5] *= -1.;

  double X1_u[4], k1_u[4];
  X1_u[0] = photon_u[0];
  X1_u[1] = photon_u[1];
  X1_u[2] = photon_u[2];
  X1_u[3] = photon_u[3];
  k1_u[0] = photon_u[4];
  k1_u[1] = photon_u[5];
  k1_u[2] = photon_u[6];
  k1_u[3] = photon_u[7];
  printf("\nNow in KS:");
  LOOP_i printf("\n%+.15e", X1_u[i]);
  LOOP_i printf("\n%+.15e", k1_u[i]);

  // Compute initial
  double k1_d[4];
  lower_index(X1_u, k1_u, k1_d);
  double BLphoton_u[8], printphoton_u[8];
  norm = four_velocity_norm(X1_u, k1_u);

// TRANSFORM BACK TO BL TO COMPUTE E, L, Q
#if (metric == KS || metric == MKS)

  KS_to_BL_u(photon_u, BLphoton_u);
  printf("\nBack to BL");
  LOOP_i printf("\n%+.15e", BLphoton_u[i]);
  LOOP_i printf("\n%+.15e", BLphoton_u[i + 4]);
  LOOP_i {
    X1_u[i] = BLphoton_u[i];
    k1_u[i] = BLphoton_u[i + 4];
  }
  BL_lower_index(X1_u, k1_u, k1_d);

#endif

  Q = k1_d[2] * k1_d[2] + cos(X1_u[2]) * cos(X1_u[2]) *
                              (a * a * (-k1_d[0] * k1_d[0]) +
                               k1_d[3] * k1_d[3] / sin(X1_u[2]) / sin(X1_u[2]));
  E = -k1_d[0];
  L = k1_d[3];

  E_init = E;
  L_init = L;
  Q_init = Q;

  printf("\nE_i = %+.15e", E);
  printf("\nL_i = %+.15e", L);
  printf("\nQ_i = %+.15e", Q);

  printf("\nNorm = %+.15e\n", norm);

  // Current r-coordinate
  r_current = logscale ? exp(photon_u[1]) : photon_u[1];

  // Reset lambda and steps
  lambda = 0.;
  steps = 0;

  // MAIN COMPUTATION
  ///////////////////

  // Trace light ray until it reaches the event horizon or the outer
  // cutoff, or steps > max_steps
#if (metric == BL || metric == MBL || metric == DM)

  // Stop condition for BL coords
  while ((r_current > cutoff_inner) && (r_current < cutoff_outer)) {
    //                   (steps < max_steps) && photon_u[0] < t_final){

#elif (metric == KS || metric == MKS)

  // Stop condition for KS coords
  while ((r_current < cutoff_outer) && steps < 1e9) {

#endif

    // Current photon position/wave vector
    LOOP_i {
      X_u[i] = photon_u[i];
      k_u[i] = photon_u[i + 4];
    }
    norm = four_velocity_norm(X_u, k_u);

    // double g_dd[4][4], g_uu[4][4];
    // metric_dd(X_u, g_dd);

    // renormalize??
    /*
    double aaa = g_dd[0][0];
    double bbb = 2. * g_dd[3][0] * k_u[3];
    double ccc = g_dd[1][1] * k_u[1] * k_u[1] + g_dd[2][2] * k_u[2]
                 * k_u[2] + g_dd[3][3] * k_u[3] * k_u[3];
    double deln=sqrt(fabs(bbb*bbb-4.*aaa*ccc));
    k_u[0]=(-bbb-deln)/2./aaa;

    double kt   = g_dd[0][0] * k_u[0] + g_dd[0][3] * k_u[3];
    double kr   = g_dd[1][1] * k_u[1];
    double kth  = g_dd[2][2] * k_u[2];
    double kphi = g_dd[3][3] * k_u[3] + g_dd[3][0] * k_u[0];
    */

    // Compute adaptive step size
    dlambda_adaptive = EPS; // stepsize(X_u, k_u);

    LOOP_i {
      printphoton_u[i] = photon_u[i];
      printphoton_u[i + 4] = photon_u[i + 4];
    }

#if (metric == KS || metric == MKS)

    KS_to_BL_u(photon_u, BLphoton_u);
    // printf("\nBack to BL");
    // LOOP_i printf("\n%+.15e", BLphoton_u[i]);
    // LOOP_i printf("\n%+.15e", BLphoton_u[i+4]);
    LOOP_i {
      printphoton_u[i] = BLphoton_u[i];
      printphoton_u[i + 4] = BLphoton_u[i + 4];
    }

#endif

    // Print current position, velocity, energy, momentum, and constants
    if (steps % 1 == 0)
      fprintf(GRfile,
              "%+.15e \t%+.15e \t%+.15e \t%+.15e \t%+.15e \t%+.15e \t%+.15e "
              "\t%+.15e \t%+.15e \t%+.15e \t%+.15e \t%+.15e \t%+.15e \n",
              lambda, printphoton_u[0],
              (logscale ? exp(printphoton_u[1]) : printphoton_u[1]),
              printphoton_u[2], printphoton_u[3], printphoton_u[4],
              printphoton_u[5], printphoton_u[6], printphoton_u[7], E, L, Q,
              norm);

      // Advance ray/particle
#if (int_method == RK4)

    rk4_step(photon_u, &f_geodesic, dlambda_adaptive);

#elif (int_method == VER)

    verlet_step(photon_u, &f_geodesic, dlambda_adaptive);

#endif

    // Put current photon into X1_u & k1_u
    X1_u[0] = photon_u[0];
    X1_u[1] = photon_u[1];
    X1_u[2] = photon_u[2];
    X1_u[3] = photon_u[3];
    k1_u[0] = photon_u[4];
    k1_u[1] = photon_u[5];
    k1_u[2] = photon_u[6];
    k1_u[3] = photon_u[7];
    //      printf("\nNORM IS %+.15e", four_velocity_norm(X1_u, k1_u));

    // Compute the constants of motion
    lower_index(X1_u, k1_u, k1_d);
    //      norm = four_velocity_norm(X1_u, k1_u);
    //      Q = k1_d[2] * k1_d[2] + cos(X1_u[2]) * cos(X1_u[2]) * (a * a *
    //          (-k1_d[0] * k1_d[0]) + k1_d[3] * k1_d[3] / sin(X1_u[2]) /
    //          sin(X1_u[2]));
    //      E = -k1_d[0];
    //      L = k1_d[3];

    // Advance (affine) parameter lambda
    lambda += fabs(dlambda_adaptive);
    r_current = logscale ? exp(photon_u[1]) : photon_u[1];

#if (metric == KS || metric == MKS)

    KS_to_BL_u(photon_u, BLphoton_u);
    t_current = BLphoton_u[0];

#endif

    steps++;
  }

  printf("\nAfter integration");
  LOOP_i printf("\n%+.15e", photon_u[i]);
  LOOP_i printf("\n%+.15e", photon_u[i + 4]);

// TEST: TRANSFORM BACK TO BL
#if (metric == KS || metric == MKS)

  KS_to_BL_u(photon_u, BLphoton_u);
  printf("\nBack to BL");
  LOOP_i printf("\n%+.15e", BLphoton_u[i]);
  LOOP_i printf("\n%+.15e", BLphoton_u[i + 4]);
  LOOP_i {
    X1_u[i] = BLphoton_u[i];
    k1_u[i] = BLphoton_u[i + 4];
  }
  BL_lower_index(X1_u, k1_u, k1_d);

#endif

  Q = k1_d[2] * k1_d[2] + cos(X1_u[2]) * cos(X1_u[2]) *
                              (a * a * (-k1_d[0] * k1_d[0]) +
                               k1_d[3] * k1_d[3] / sin(X1_u[2]) / sin(X1_u[2]));
  E = -k1_d[0];
  L = k1_d[3];

  printf("\nE_f = %+.15e", E);
  printf("\nL_f = %+.15e", L);
  printf("\nQ_f = %+.15e", Q);

  printf("\nL1 norm:");
  printf("\nE error is: %+.15e ", fabs(E - E_init) / E_init);
  printf("\nL error is: %+.15e ", fabs(L - L_init) / L_init);
  printf("\nQ error is: %+.15e ", fabs(Q - Q_init) / Q_init);

  // CLOSE & EXIT
  ///////////////

  fclose(GRfile);

  return 0;
}
