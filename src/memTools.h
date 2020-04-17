/** \file
 *
 * \author hhh
 * \date Fri 27 Mar 2020 02:38:45 PM CET
 */

#ifndef MEMTOOLS_H
#define MEMTOOLS_H

#include "cgnslib.h"

long **dyn2DintArray(long I, long J);
cgsize_t **dyn2DcgsizeArray(long I, long J);
double **dyn2DdblArray(long I, long J);
long ***dyn3DintArray(long I, long J, long K);
double ***dyn3DdblArray(long I, long J, long K);
double ****dyn4DdblArray(long I, long J, long K, long L);
char **dynStringArray(long I, long J);

#endif
