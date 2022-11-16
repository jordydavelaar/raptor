/*
 * model_global_vars.h
 *
 * Written by J. Davelaar 2022
 *
 */

#ifndef MODEL_GLOBAL_VARS_H
#define MODEL_GLOBAL_VARS_H

extern double ****p;

extern double L_unit, T_unit;
extern double RHO_unit, U_unit, B_unit;
extern double Ne_unit, Thetae_unit;

extern double R0, a, Q, hslope;
extern double stopx[4], startx[4], dx[4];

extern double ***Xgrid;
extern struct block *block_info;
#endif // MODEL_GLOBAL_VARS_H
