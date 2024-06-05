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

auto jnu_model = Model::load("jnu-v5.model");

void int_jnu_ml(){
    // Initialize model.
    auto jnu_model = Model::load("jnu-v5.model");
}

double jnu_ml(double nu_p,struct GRMHD  modvar, double pitch_ang) {
    double nuc = modvar.B * ELECTRON_CHARGE/(2*M_PI*ELECTRON_MASS*SPEED_OF_LIGHT);

    double mu=cos(pitch_ang);
    double X[3]={log10(nu_p/nuc),log10(modvar.theta_e),mu};

    double xmin[3] = {2.7369195028346725e-05 , -1.9999985011285943,-0.9999998523544762}; 
    double xmax[3] = {8.99996224216519,3.999985633637197,1.9999992996211267};
    double ymin = -298.9876181337542;
    double ymax = 298.916311406278;

    X[0]=(X[0]-xmin[0])/xmax[0];
    X[1]=(X[1]-xmin[1])/xmax[1];
    X[2]=(X[2]-xmin[2])/xmax[2];

    // Create a 1D Tensor on length 10 for input data.
    Tensor in{3};
    in.data_ = {X[0], X[1], X[2]};

    if(X[0]<0 || X[0]>1){
	return 0;//j_I_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);;
	}
    if(X[1]<0 || X[1]>1){ 
	return 0;//j_I_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);;
    }
    if(X[2]<0 || X[2]>1){ 
	return 0;//j_I_thermal(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);;
    }

    // Run prediction.
    Tensor out = jnu_model(in);
    double jfac = ELECTRON_CHARGE*ELECTRON_CHARGE * nuc /SPEED_OF_LIGHT;
    double jnu = modvar.n_e*pow(10.,out.data_[0]*ymax+ymin)*jfac;

    return jnu;
}
