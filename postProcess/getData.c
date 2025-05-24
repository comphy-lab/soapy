/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025
*/

#include "utils.h"
#include "output.h"

scalar f[];
vector u[];

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

scalar T[];
scalar * list = NULL;

/**
   * @brief Entry point for simulation snapshot extraction and processing.
   *
   * This function validates command-line arguments and orchestrates the simulation
   * data restoration, derivative and velocity computations, and interpolation onto
   * a grid. It expects the program name followed by exactly six parameters:
   * a filename, the lower bounds (xmin and ymin), the upper bounds (xmax and ymax), and
   * the number of divisions along the y-axis (ny). If the argument count is incorrect,
   * an error message and usage instructions are printed to stderr and the program exits
   * with a status of 1.
   *
   * @param a The total number of command-line arguments.
   * @param arguments Array of command-line argument strings, where arguments[0] is the
   * program name and the remaining elements provide the required simulation parameters.
   *
   * @return int Exit status of the program (1 on error, 0 on success).
   */
  int main(int a, char const *arguments[])
{
  if (a != 7) {
    fprintf(stderr, "Error: Expected 6 arguments\n");
    fprintf(stderr, "Usage: %s <filename> <xmin> <ymin> <xmax> <ymax> <ny>\n", arguments[0]);
    return 1;
  }

  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);

  list = list_add (list, T);

  /*
  Actual run and codes!
  */
  restore (file = filename);

  double maxT = statsf(T).max;
  foreach(){
    T[] *= 1e0/maxT;
  }

  FILE * fp = ferr;
  Deltay = (double)((ymax-ymin)/(ny));
  // fprintf(ferr, "%g\n", Deltay);
  nx = (int)((xmax - xmin)/Deltay);
  // fprintf(ferr, "%d\n", nx);
  Deltax = (double)((xmax-xmin)/(nx));
  // fprintf(ferr, "%g\n", Deltax);
  len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
}
