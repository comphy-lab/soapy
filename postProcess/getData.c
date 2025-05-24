/**
# Simulation Data Extraction

Extract and interpolate field data from Basilisk snapshots onto a regular grid.

- **Author**: Vatsal Sanjay (vatsalsy@comphy-lab.org)
- **Lab**: CoMPhy Lab, Physics of Fluids Department
- **Updated**: March 8, 2025

## Usage
```bash
./getData <filename> <xmin> <ymin> <xmax> <ymax> <ny>
```

Output: Grid data (x, y, field_values) to stderr
*/

#include "utils.h"
#include "output.h"

/** VOF and velocity fields */
scalar f[];
vector u[];

/** Command-line parameters */
char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay;

/** Temperature field and scalar list */
scalar T[];
scalar * list = NULL;

/**
## Main Function

Extract simulation data and interpolate onto a regular grid.

### Parameters
- `filename`: Simulation snapshot file
- `xmin`, `ymin`: Lower bounds of extraction region
- `xmax`, `ymax`: Upper bounds of extraction region  
- `ny`: Number of grid points in y-direction

### Process
1. Parse arguments and validate input
2. Load simulation and normalize T field
3. Create uniform grid with spacing Deltay
4. Interpolate fields onto grid points
5. Output grid data to stderr
*/
int main(int a, char const *arguments[]) {
  if (a != 7) {
    fprintf(stderr, "Error: Expected 6 arguments\n");
    fprintf(stderr, "Usage: %s <filename> <xmin> <ymin> <xmax> <ymax> <ny>\n", 
            arguments[0]);
    return 1;
  }

  sprintf(filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);

  list = list_add(list, T);

  /* Load simulation */
  restore(file = filename);

  /* Normalize T by maximum value */
  double maxT = statsf(T).max;
  foreach() {
    T[] *= 1e0/maxT;
  }

  FILE * fp = ferr;
  
  /* Calculate grid parameters */
  Deltay = (double)((ymax-ymin)/(ny));
  nx = (int)((xmax - xmin)/Deltay);
  Deltax = (double)((xmax-xmin)/(nx));
  len = list_len(list);
  
  /* Allocate and fill interpolation matrix */
  double ** field = (double **) matrix_new(nx, ny+1, len*sizeof(double));
  
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list) {
        field[i][len*j + k++] = interpolate(s, x, y);
      }
    }
  }

  /* Output interpolated data */
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf(fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list) {
        fprintf(fp, " %g", field[i][len*j + k++]);
      }
      fputc('\n', fp);
    }
  }
  
  fflush(fp);
  matrix_free(field);
}