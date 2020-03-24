/*
 * boundary.h
 *
 * Created: Tue 24 Mar 2020 10:02:10 AM CET
 * Author : hhh
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <stdbool.h>

#include "main.h"

typedef struct boundary_t boundary_t;
struct boundary_t {
	int BCtype;	/* Boundary type:
			 * 1: Slip Wall
			 * 2: Navier-Stokes Wall
			 * 3: Supersonic Inflow
			 * 4: Supersonic Outflow
			 * 5: Characteristic
			 * 6: Exact Solution */
	int BCid;

	int exactFunc;

	double pVar[NVAR];

	double wallVelocity;
	bool isAdiabatic;
	bool isTemperaturePrescribed;	/* true: wall temperature prescribed
					 * false: heat flux prescribed */
	double temperature;
	double heatFlux;

	double *connection;

	boundary_t *next;
};
extern boundary_t *firstBC;
extern int nBC;
extern bool isPeriodic;

void initBoundary(void);

#endif
