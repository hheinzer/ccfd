/*
 * equation.h
 *
 * Created: Tue 24 Mar 2020 08:23:21 AM CET
 * Author : hhh
 */

#ifndef EQUATION_H
#define EQUATION_H

#include <stdbool.h>

extern double pi;

extern bool doCalcSource;
extern double R;	/* gas constant */
extern double gam;	/* ratio of specific heat */
extern double gam1;	/* gam - 1 */
extern double gam2;	/* gam - 2 */
extern double gam1q;	/* 1 / (gam - 1) */
extern double cp;	/* specific heat capacity */
extern double Pr;	/* Prandtl number */
extern double mu;	/* viscosity */

extern int iFlux;

extern int intExactFunc;
extern int sourceFunc;
extern double sqrt3q, sqrt3, sqrt2;

void initEquation(void);

#endif
