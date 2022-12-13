/*
 * model_functions.h
 *
 * harmpi model for RAPTOR
 *
 * Based on public harmpi code: https://github.com/atchekho/harmpi
 *
 * Authors Anna Chashkina, Jordy Davelaar
 *
 */

#ifndef MODEL_FUNCTIONS_H
#define MODEL_FUNCTIONS_H

double interp_scalar(double ***var, int i, int j, int k, double coeff[4]);

#endif // RAPTOR_HARM_MODEL_H
