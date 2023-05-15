/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f2[];
char filename[80];
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  #if TREE
    f2.prolongation = fraction_refine;
  #endif
  boundary((scalar *){f2});
  FILE * fp = ferr;
  output_facets(f2,fp);
  fflush (fp);
  fclose (fp);
}
