/**
# Facet Extraction Utility

Extract interface facets from Basilisk simulation files.

- **Author**: Vatsal Sanjay (vatsalsy@comphy-lab.org)
- **Lab**: CoMPhy Lab, Physics of Fluids Department
- **Updated**: March 8, 2025

## Usage
```bash
./getFacets <simulation_file>
```

Output: Facet data to stderr
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

/** VOF fraction field */
scalar f[];

/** Input filename buffer */
char filename[80];

/**
## Main Function

Extract facets from a saved simulation state.

### Parameters
- `a`: Argument count
- `arguments`: Command-line arguments
  - `arguments[1]`: Simulation file path

### Process
1. Load simulation file
2. Extract facets from VOF field
3. Output to stderr
*/
int main(int a, char const *arguments[]) {
  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);
  
  FILE * fp = ferr;
  output_facets(f, fp);
  
  fflush(fp);
  fclose(fp);
}