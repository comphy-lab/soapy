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