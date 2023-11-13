#include <math.h>
#include <stdio.h>
#include "mfobj.h"

#define SQUARE(x) ((x) * (x))

//-----------------------------------------------------------------------------
// MFObj function
// - n is the dimension of the data
// - point is the location where the function will be evaluated
// - arg contains the parameters of the function
// More details on the function at http://www.sfu.ca/%7Essurjano/ackley.html
//-----------------------------------------------------------------------------

void mfobj_fun(int n, point_t *point, const void *arg) {

  // cast the void pointer to what we expect to find
  const mfobj_param_t *mf_params = (const mfobj_param_t *)arg;
  int m = mf_params->ysquare_size;
  // cost function computation
  double result[m]; // Array of the size of y
  double obj = 0;
  for (int j = 0; j < m; j++) {
      result[j] = 0;
      for (int i = 0; i < n-1; i++) {
          result[j] = result[j] + point->x[i] * mf_params->Asmall[j][mf_params->dicPos[i]];
      }
    //printf("Ysquare: %f\n", mf_params->ysquare[j]);
    //printf("Result: %f\n", sqrt(SQUARE(result[j]) + SQUARE(point->x[n-1])));
      obj = obj + SQUARE(mf_params->ysquare[j] - sqrt(SQUARE(result[j]) + SQUARE(point->x[n-1])));
  }

  //Penalize negative x values on the objective function
    double totWeight = 0.0;
    for (int i = 0; i < n; i++) {
        if (point->x[i] < 0) {
            obj = obj - log(point->x[i]);
        }
        if (i<n-1){
            totWeight = point->x[i] + totWeight;
        }

    }

    //Penalize if the sum of the weights is not 1
    obj = obj;// + 2*SQUARE(totWeight - 1);

  // final result
  point->fx = obj;
}
