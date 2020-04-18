/** \file
 *
 * \author hhh
 * \date Fri 27 Mar 2020 02:28:25 PM CET
 */

#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

extern int icType;
extern int nDomains;
extern int *domainID;
extern double rp1Dinterface;
extern double alpha;
extern double **refState;

void initInitialCondition(void);
void setInitialCondition(void);
void freeInitialCondition(void);

#endif
