/** \file
 *
 * \author hhh
 * \date Sat 28 Mar 2020 09:45:50 PM CET
 */

#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "main.h"

void primCons(const double pVar[NVAR], double cVar[NVAR]);
void consPrim(const double cVar[NVAR], double pVar[NVAR]);
void consChar(double cVar[NVAR], double charac[3], double pVarRef[NVAR]);
void charCons(double charac[3], double cVar[NVAR], double pVarRef[NVAR]);

#endif
