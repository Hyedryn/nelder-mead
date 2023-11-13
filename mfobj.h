#ifndef MFOBJ_H
#define MFOBJ_H

#include "point.h"
#include <stddef.h>

//-----------------------------------------------------------------------------
// Implementation of a cost function f : R^n->R compatible with the fun_t  
// interface defined in cost.h; here we use the Ackley Function as it allows 
// us to demonstrate the use of the optional fixed arguments.
//-----------------------------------------------------------------------------

typedef struct {
  int ysquare_size;
  double* ysquare;
  double** Asmall;
  int* dicPos;
} mfobj_param_t;

void mfobj_fun(int, point_t *, const void *);

#endif // MFOBJ_H
