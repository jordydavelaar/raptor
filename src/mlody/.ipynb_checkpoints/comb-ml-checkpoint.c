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

auto model  = Model::load("comb_final_compact.model");

double xmin[3];
double xmax[3];

double ymin[8];
double ymax[8];

void init_ml_model(){

    FILE *input;

    input = fopen("minmax.txt", "r");
    if (input == NULL) {
        fprintf(stderr, "Can't read file %s! Aborting", "minmax.txt");
        exit(1);
    }

    // Model parameters
    for(int i=0;i<3;i++){
	fscanf(input, "%lf %lf", &xmin[i],&xmax[i]);
	 
   }

    for(int i=0;i<8;i++){
	fscanf(input, "%lf %lf", &ymin[i],&ymax[i]);
    }


}

void comb_ml(double nu_p,struct GRMHD  modvar, double pitch_ang, 
             double *jI, double *aI, double *jQ, double *aQ, double *rQ, 
             double *jV, double *aV, double *rV){
    double nuc = modvar.B * ELECTRON_CHARGE/(2*M_PI*ELECTRON_MASS*SPEED_OF_LIGHT);

    double mu=cos(pitch_ang);
    double X[3]={log10(nu_p/nuc),log10(modvar.theta_e),mu};

    X[0]=(X[0]-xmin[0])/xmax[0];
    X[1]=(X[1]-xmin[1])/xmax[1];
    X[2]=(X[2]-xmin[2])/xmax[2];


    // Create a 1D Tensor on length 10 for input data.
    Tensor in{3};
    in.data_ = {X[0], X[1], X[2]};

    if(X[0]<-1 || X[0]>1 || X[1]<-1 || X[1]>1 || X[2]<-1 || X[2]>1){ 
	    *jI = j_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
	    *jQ = j_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
	    *jV = j_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

	    *rQ = rho_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);
	    *rV = rho_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang);

	    *aI = a_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jI);
	    *aQ = a_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jQ);
	    *aV = a_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, *jV);

	return;
    }

    // Run prediction.
    Tensor out = model(in);
    double jfac = ELECTRON_CHARGE*ELECTRON_CHARGE * nuc /SPEED_OF_LIGHT;
    double afac = ELECTRON_CHARGE*ELECTRON_CHARGE/(nu_p*ELECTRON_MASS*SPEED_OF_LIGHT);

    *jI = modvar.n_e*pow(10.,out.data_[0]*ymax[0]+ymin[0])*jfac;
    *aI = modvar.n_e*pow(10.,out.data_[1]*ymax[1]+ymin[1])*afac;

    *jQ = modvar.n_e*pow(10.,out.data_[2]*ymax[2]+ymin[2])*jfac;
    *aQ = modvar.n_e*pow(10.,out.data_[3]*ymax[3]+ymin[3])*afac;
    *rQ = -modvar.n_e *sign(out.data_[4])*powf(out.data_[4]*ymax[4],4.0);

    *jV = mu*modvar.n_e*pow(10.,out.data_[5]*ymax[5]+ymin[5])*jfac;
    *aV = mu*modvar.n_e*pow(10.,out.data_[6]*ymax[6]+ymin[6])*afac;
    *rV = modvar.n_e *sign(out.data_[7])*powf(out.data_[7]*ymax[7],4.0);
// 	fprintf(stderr,"%e %e %e %e %e %e %e %e %e\n",*jI,*aI,*jQ,*aQ,*rQ,*jV,*aV,*rV);
}


