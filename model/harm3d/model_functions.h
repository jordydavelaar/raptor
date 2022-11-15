/*
 * model_functions.h
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
 */

#ifndef MODEL_FUNCTIONS_H
#define MODEL_FUNCTIONS_H

double interp_scalar(double ***var, int i, int j, int k, double coeff[4]);

#endif // RAPTOR_HARM_MODEL_H
