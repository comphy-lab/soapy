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
- `argc`: Argument count
- `argv`: Command-line arguments
  - `argv[1]`: Simulation file path

### Process
1. Validate arguments
2. Load simulation file
3. Extract facets from VOF field
4. Output to stderr
*/
int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <simulation_file>\n", argv[0]);
    return 1;
  }
  
  snprintf(filename, sizeof(filename), "%s", argv[1]);
  restore(file = filename);
  
  FILE * fp = ferr;
  output_facets(f, fp);
  
  fflush(fp);
  fclose(fp);
}