/*
 * analyze.h
 *
 * Created: Sun 29 Mar 2020 05:51:23 PM CEST
 * Author : hhh
 */

#ifndef ANALYZE_H
#define ANALYZE_H

#include <stdlib.h>
#include <stdbool.h>

#include "main.h"
#include "boundary.h"
#include "mesh.h"

typedef struct times_t times_t;
struct times_t {
	double overall;
	double io;
	double flux;
	double spaceRec;
	double ck;
	double timeInt;
};

typedef struct wing_t wing_t;
struct wing_t {
	double refLength;
	int wallId;
	double cl, cd;
	boundary_t *wingBC;
	sidePtr_t *firstPressureSide;
	sidePtr_t *firstSuctionSide;
};

typedef struct recordPoint_t recordPoint_t;
struct recordPoint_t {
	int nPoints;
	double **x;
	elem_t **elem;
	FILE **ioFile;
};

extern bool doCalcWing;
extern wing_t wing;
extern recordPoint_t recordPoint;
extern bool hasExactSolution;

void initAnalyze(void);
void analyze(double time, long iter, double *resIter);
void calcErrors(double time);
void globalResidual(double dt, double resIter[NVAR + 2]);

#endif
