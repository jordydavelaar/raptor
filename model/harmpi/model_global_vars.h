/*
 * model_global_vars.h
 *
 * harmpi model for RAPTOR
 *
 * Based on public harmpi code: https://github.com/atchekho/harmpi
 *
 * Authors Anna Chashkina, Jordy Davelaar
 *
 */

#ifndef MODEL_GLOBAL_VARS_H
#define MODEL_GLOBAL_VARS_H

extern double ****p;

extern int N1, N2, N3;

extern double R_HIGH, R_LOW, gam;

extern double R0, Rin, Rout, a, hslope;
extern double startx[NDIM], stopx[NDIM], dx[NDIM];

extern double L_unit, T_unit;
extern double RHO_unit, U_unit, B_unit;
extern double Ne_unit, Thetae_unit;

#endif // MODEL_GLOBAL_VARS_H
