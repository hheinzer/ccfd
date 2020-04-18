/** \file
 *
 * \brief Contains the structure definitions of `wing_t` and `recordPoint_t`
 *
 * \author hhh
 * \date Sun 29 Mar 2020 05:51:23 PM CEST
 */

#ifndef ANALYZE_H
#define ANALYZE_H

typedef struct wing_t wing_t;
typedef struct recordPoint_t recordPoint_t;

#include <stdlib.h>
#include <stdbool.h>

#include "main.h"
#include "boundary.h"
#include "mesh.h"

/**
 * \brief Collection of all necessary values for the calculation of CL and CD
 */
struct wing_t {
	double refLength;		/**< reference length */
	int wallId;			/**< BCid of the wall that represents the wing */
	double cl;			/**< lift coefficient */
	double cd;			/**< drag coefficient */
	boundary_t *wingBC;		/**< pointer to the BC of the wing */
	sidePtr_t *firstPressureSide;	/**< pointer to the first pressure side */
	sidePtr_t *firstSuctionSide;	/**< pointer to the first suction side */
};

/**
 * \brief Recording point structure, used to output flow field at specific points
 */
struct recordPoint_t {
	int nPoints;			/**< number of recording points */
	double **x;			/**< `nPoints`x`NDIM` array of RP coordinates */
	elem_t **elem;			/**< array of elements that contain a RP */
	FILE **ioFile;			/**< array of output file pointers */
};

extern bool doCalcWing;
extern wing_t wing;
extern recordPoint_t recordPoint;
extern bool hasExactSolution;

void initAnalyze(void);
void analyze(double time, long iter, double *resIter);
void calcErrors(double time);
void globalResidual(double resIter[NVAR + 2]);
void freeAnalyze(void);

#endif
