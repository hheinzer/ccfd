/** \file
 *
 * \author hhh
 * \date Tue 24 Mar 2020 08:23:21 AM CET
 */

#ifndef EQUATION_H
#define EQUATION_H

#include <stdbool.h>

extern double pi;

extern bool doCalcSource;
extern double R;
extern double gam;
extern double gam1;
extern double gam2;
extern double gam1q;
extern double cp;
extern double Pr;
extern double mu;

extern int iFlux;

extern int intExactFunc;
extern int sourceFunc;
extern double sqrt2;
extern double sqrt3;
extern double sqrt3q;

void initEquation(void);

#endif
