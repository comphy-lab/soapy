/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsy@comphy-lab.org
# CoMPhy Lab
# Physics of Fluids Department
# Last updated: Mar 8, 2025
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];
char filename[80];

int main(int a, char const *arguments[]){
  sprintf(filename, "%s", arguments[1]);

  restore (file = filename);

  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
