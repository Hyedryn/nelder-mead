#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mfobj.h"
#include "nelder_mead.h"

#define PI 3.1415926535897932384626433832795

//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------

int main(int argc, const char *argv[]) {
  if (argc == 1) {
    printf("%s: error: not enough inputs \n", argv[0]);
    return 0;
  }

  // reading initial point from command line
  // 0.5 0.5 1 380 512
  const int n = argc - 3;
  point_t start;
  start.x = malloc(n * sizeof(double));
  for (int i = 0; i < n-1; i++) {
    start.x[i] = atof(argv[i + 1]);
  }
  int num_atoms = atoi(argv[n+1]);
  int num_y = atoi(argv[n+2]);
    printf("Num atoms: %d\n", num_atoms);
    printf("Num y: %d\n", num_y);

  // optimisation settings
  optimset_t optimset;
  optimset.tolx = 0.001;
  optimset.tolf = 0.001;
  optimset.max_iter = 1000;
  optimset.max_eval = 1000;
  optimset.verbose = 0;

  // cost function parameters
  mfobj_param_t mf_params;
  mf_params.ysquare_size = num_y;
  mf_params.dicPos = malloc(sizeof(int) * (n-1));
  mf_params.dicPos[0] = 0;
  mf_params.dicPos[1] = num_atoms;
  mf_params.dicPos[2] = num_atoms + num_atoms - 25;


// Allocate memory for mf_params.Asmall
   double **Asmall = calloc(num_y, sizeof(double *));
   if (Asmall == NULL) {
        fprintf(stderr, "Memory allocation error for rows\n");
        return 1;
   }

   for (int i = 0; i < num_y; i++) {
        Asmall[i] = calloc(( (2*num_atoms) +1), sizeof(double));
        if (Asmall[i] == NULL) {
            fprintf(stderr, "Memory allocation error for columns\n");
            return 1;
        }
    }


    FILE *fileAsmall = fopen("dicA.txt", "r");
    if (fileAsmall == NULL) {
        perror("Error opening file");
        return 1;
    }

    for (int i = 0; i < num_y; i++) {
        for (int j = 0; j < (2*num_atoms)+1; j++) {
            if (fscanf(fileAsmall, "%lf", &Asmall[i][j]) != 1) {
                fprintf(stderr, "Error reading from file");
                return 1;
            }
        }
    }

    double *ysquare = malloc(num_y * sizeof(double));
    if (ysquare == NULL) {
        fprintf(stderr, "Memory allocation error for ysquare\n");
        return 1;
    }
    FILE *filey = fopen("y.txt", "r");
    if (filey == NULL) {
        perror("Error opening file");
        return 1;
    }
    for (int i = 0; i < num_y; i++) {
        if (fscanf(filey, "%lf", &ysquare[i]) != 1) {
            fprintf(stderr, "Error reading from file");
            return 1;
        }
    }

    // Put y to square
    for (int i=0; i<num_y; i++) {
        ysquare[i] = ysquare[i]/3300.0;
        ysquare[i] = ysquare[i]*ysquare[i];
    }

    mf_params.ysquare = ysquare;
    mf_params.Asmall = Asmall;


    //Print an element of Asmall
    printf("Asmall[0][0]: %f\n", Asmall[0][0]);
    //Print an element of mf_params.Asmall
    printf("mf_params.Asmall[0][0]: %f\n", mf_params.Asmall[0][0]);

    //Print an element of ysquare
    printf("ysquare[0]: %f\n", ysquare[0]);
    //Print an element of mf_params.ysquare
    printf("mf_params.ysquare[0]: %f\n", mf_params.ysquare[0]);



  // sum of elem of ysquare
  double init_cost=0;
  for (int k=0; k< mf_params.ysquare_size; k++) {
      init_cost = init_cost + mf_params.ysquare[k];
  }

  printf("Init cost: %f\n", init_cost);

  double best_cost = init_cost;
  int best_i=-1;
  int best_j=-1;
  point_t best_solution;
  clock_t startTime = clock();
  for (int i = 0; i < num_atoms; i++) {
    for (int j = 0; j < num_atoms; j++) {
      // call optimization method
      mf_params.dicPos[0] = i;
      mf_params.dicPos[1] = num_atoms + j;
      point_t solution;
      nelder_mead(n, &start, &solution, &mfobj_fun, &mf_params, &optimset);

      // Print solution.fx
      //printf("Solution.fx: %f , j: %d\n", solution.fx, j);
      if ((solution.fx < best_cost) || (i==0 && j==0)) {
        best_cost = solution.fx;
        best_solution = solution;
        best_i = i;
        best_j = j;
      }
    }
      if (i%5 == 0) {
          clock_t middleTime = clock();
          printf("Best solution so far: %f\n", best_cost);
          printf("Best i and j so far: %d %d\n", best_i, best_j);
          printf("Solution x so far: %f %f %f %f\n", best_solution.x[0], best_solution.x[1], best_solution.x[2], best_solution.x[3]);
          double cpu_time_used = ((double) (middleTime - startTime)) / CLOCKS_PER_SEC;
          double iter = i*num_atoms;
          double iterBySecond = iter/cpu_time_used;
          double eta = (num_atoms*num_atoms - iter)/iterBySecond;
          printf("Iteration: %d (% 0.2f iter/s) ETA: % 0.2f \n", i, iterBySecond, eta);

      }
  }
  clock_t endTime = clock();
    double cpu_time_used = ((double) (endTime - startTime)) / CLOCKS_PER_SEC;
    printf("Total duration: %f\n", cpu_time_used);



  // evaluate and print starting point
  printf("Initial point\n");
  mfobj_fun(n, &start, &mf_params);
  print_point(n, &start);
  // print solution
  printf("Solution\n");
  printf("I and J selected: %d %d\n", best_i, best_j);
  print_point(n, &best_solution);

    FILE *fileW = fopen("w.txt", "w");
    if (fileW == NULL) {
        perror("Error opening file");
        return -1;
    }
    for (int i = 0; i < n; i++) {
        fprintf(fileW, "%f ", best_solution.x[i]);
    }
    fclose(fileW);

    FILE *fileInd = fopen("ind.txt", "w");
    if (fileInd == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(fileInd, "%d %d 0", best_i, best_j);
    fclose(fileInd);

    FILE *fileObj = fopen("obj.txt", "w");
    if (fileObj == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(fileObj, "%f", best_solution.fx);
    fclose(fileObj);

  // free memory
  free(start.x);
  free(best_solution.x);

  return 0;
}
