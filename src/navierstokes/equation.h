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
extern double gamma;	/* ratio of specific heat */
extern double gamma1;	/* gamma - 1 */
extern double gamma2;	/* gamma - 2 */
extern double gamma1q;	/* 1 / (gamma - 1) */
//extern double CFL;	/* CFL number */
extern double cp;	/* specific heat capacity */
extern double Pr;	/* Prandtl number */
extern double ReRef;	/* Reference Reynolds number */
extern double mu;	/* viscosity */

extern int iFlux;

extern int intExactFunc;
extern int sourceFunc;
extern double sqrt3q, sqrt3, sqrt2;

void initEquation(void);

#endif
