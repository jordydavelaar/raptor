#include "model-ml.h"

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <iostream>
#include <math.h>

using keras2cpp::Model;
using keras2cpp::Tensor;

auto anu_model = Model::load("anu-v5.model");


double anu_ml(double nu_p,struct GRMHD  modvar, double pitch_ang) {
    double nuc = modvar.B * ELECTRON_CHARGE/(2*M_PI*ELECTRON_MASS*SPEED_OF_LIGHT);

    double mu=cos(pitch_ang);
    double X[3]={log10(nu_p/nuc),log10(modvar.theta_e),mu};

    double xmin[3] = {2.7369195028346725e-05, -1.9999985011285943, -0.9999998523544762};
    double xmax[3] = {8.99996224216519, 3.999985633637197, 1.9999992996211267};
    double ymin = -298.99997021704286;
    double ymax = 300.35167346230116;

    X[0]=(X[0]-xmin[0])/xmax[0];
    X[1]=(X[1]-xmin[1])/xmax[1];
    X[2]=(X[2]-xmin[2])/xmax[2];

    // Create a 1D Tensor on length 10 for input data.
    Tensor in{3};
    in.data_ = {X[0], X[1], X[2]};

    if(X[0]<0 || X[0]>1)
	return 0;
    if(X[1]<0 || X[1]>1)
	return 0;
    if(X[2]<0 || X[2]>1)
	return 0;


    // Run prediction.
    Tensor out = anu_model(in);
    double afac = ELECTRON_CHARGE*ELECTRON_CHARGE/(nu_p*ELECTRON_MASS*SPEED_OF_LIGHT);
    double anu = modvar.n_e*pow(10.,out.data_[0]*ymax+ymin)*afac;

    return anu;
}
