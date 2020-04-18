/** \file
 *
 * \author hhh
 * \date Sat 28 Mar 2020 10:16:16 AM CET
 */

#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

extern int limiter;
extern double venk_k;
extern int nVarGrad;

void spatialReconstruction(double time);

#endif
