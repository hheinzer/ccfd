/** \file
 *
 * \author hhh
 * \date Fri 27 Mar 2020 05:09:57 PM CET
 */

#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

extern int spatialOrder;
extern int fluxFunction;

void initFV(void);
void fvTimeDerivative(double time);

#endif
