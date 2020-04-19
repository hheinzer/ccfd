/** \file
 *
 * \brief Contains the structure definition of a boundary
 *
 * \author hhh
 * \date Tue 24 Mar 2020 10:02:10 AM CET
 */

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <stdbool.h>

typedef struct boundary_t boundary_t;

#include "main.h"
#include "mesh.h"

/**
 * \brief Structure that holds the information of a boundary condition
 */
struct boundary_t {
	int BCtype;			/**< boundary type:
						- 1: slip wall
						- 2: Navier-Stokes wall
						- 3: supersonic inflow
						- 4: supersonic outflow
						- 5: characteristic
						- 6: exact solution
						- 7: periodic boundary condition
						- 8: pressure outflow */
	int BCid;			/**< boundary condition ID */

	int exactFunc;			/**< exact boundary function identifier */

	double pVar[NVAR];		/**< inflow state */

	bool isAdiabatic;		/**< adiabatic wall flag */
	bool isTemperaturePrescribed;	/**< is the temperature prescribed flag */
	double temperature;		/**< wall temperature */
	double heatFlux;		/**< wall heat flux */

	double *connection;		/**< connection coordinates for periodic BC */

	boundary_t *next;		/**< pointer to next boundary condition */
};

extern boundary_t *firstBC;
extern int nBC;
extern bool isPeriodic;

void initBoundary(void);
void setBCatSides(double time);
void setBCatBarys(double time);
void boundary(side_t *aSide, double time, double int_pVar[NVAR],
		double ghost_pVar[NVAR], double x[NDIM]);
void freeBoundary(void);

#endif
