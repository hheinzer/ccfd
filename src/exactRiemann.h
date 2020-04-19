/** \file
 *
 * \date Sun 29 Mar 2020 12:35:26 PM CEST
 * \author hhh
 */

#ifndef EXACTRIEMANN_H
#define EXACTRIEMANN_H

void exactRiemann(double rhol, double rhor, double *rho,
		  double ul,   double ur,   double *u,
		  double pl,   double pr,   double *p,
		  double al,   double ar,   double s);

#endif
