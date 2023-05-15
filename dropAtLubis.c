/** Title: Drop spreading at (miscible) LUBIS
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Updated: Mon 15th May 2023
*/

// 1 is oil layer, 2 is gas bubble and 3 is water
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"

#define MINlevel 3                                              // maximum level

#define tsnap (1e-3)

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-4)                                 // error tolerance in KAPPA
#define VelErr (1e-3)                            // error tolerances in velocity
#define OmegaErr (1e-2)                            // error tolerances in velocity

// boundary conditions
f1[left] = dirichlet(1.0);
f2[left] = dirichlet(0.0);
u.t[left] = dirichlet(0.0);


double Ohf, Ohd, Ohs, hf, tmax, Ldomain, delta;
int MAXlevel;

int main(int argc, char const *argv[]) {

  if (argc < 6) {
    fprintf(ferr, "Need %d more argument(s): Oho, Ohw, Oha, hf, tmax, Ldomain, delta, MAXlevel\n", 6-argc);
    return 1;
  }
  
  Ohf = atof(argv[1]); // film Ohnesorge number
  Ohd = atof(argv[2]); // drop Ohnesorge number
  Ohs = 1e-5; // surrounding Ohnesorge number
  hf = atof(argv[3]);
  tmax = atof(argv[4]);
  Ldomain = atof(argv[5]);
  delta = 1e-2;
  MAXlevel = atoi(argv[6]);

  L0=Ldomain;
  X0=-hf*1.001; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1e0; mu1 = Ohf;
  rho2 = 1e0; mu2 = Ohd;
  rho3 = 1e-3; mu3 = Ohs;

  f1.sigma = 1.0;
  f2.sigma = 0.0; // drop and film has zero surface tension

  fprintf(ferr, "Level %d tmax %g. Ohf %3.2e, Ohd %3.2e, Ohs %3.2e, hf %3.2f\n", MAXlevel, tmax, Ohf, Ohd, Ohs, hf);
  run();
}

int refRegion(double x, double y, double z){

  return (x < 1.1*hf && y < 1.0 ? MAXlevel:
          x < 5*hf && y < 1.5 ? MAXlevel-1:
          x < 1.0 && y < 2.0 ? MAXlevel-2:
          x < 3.0 && y < 2.0 ? MAXlevel-3:
          MAXlevel-4
        );
}

event init(t = 0){
  if(!restore (file = "dump")){
    char filename1[60], filename2[60];
    /**
    Initialization for f1 and f2
    */
    sprintf(filename1, "f1_init-%3.2f.dat", delta);
    sprintf(filename2, "f2_init-%3.2f.dat", delta);

    FILE * fp1 = fopen(filename1,"rb");
    if (fp1 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename1);
      return 1;
    }
    FILE * fp2 = fopen(filename2,"rb");
    if (fp2 == NULL){
      fprintf(ferr, "There is no file named %s\n", filename2);
      return 1;
    }
    coord* InitialShape1;
    coord* InitialShape2;
    InitialShape1 = input_xy(fp1);
    fclose (fp1);
    InitialShape2 = input_xy(fp2);
    fclose (fp2);
    scalar d1[], d2[];
    distance (d1, InitialShape1);
    distance (d2, InitialShape2);
    while (adapt_wavelet_limited ((scalar *){f1, f2, d1, d2}, (double[]){1e-8, 1e-8, 1e-8, 1e-8}, refRegion).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi1[], phi2[];
    foreach_vertex(){
      phi1[] = -(d1[] + d1[-1] + d1[0,-1] + d1[-1,-1])/4.;
      phi2[] = -(d2[] + d2[-1] + d2[0,-1] + d2[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fractions in the domain. */
    fractions (phi1, f1);
    fractions (phi2, f2);
  }
  dump (file = "dump");
  // return 1;
}

scalar KAPPA1[], KAPPA2[], omega[];
event adapt(i++) {
  vorticity (u, omega);
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  foreach(){
    omega[] *= f1[]*(1-f2[]);
  }
  adapt_wavelet_limited ((scalar *){f1, f2, u.x, u.y, KAPPA1, KAPPA2, omega},
    (double[]){fErr, fErr, VelErr, VelErr, KErr, KErr, OmegaErr},
    refRegion, MINlevel);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f1[],f2[]);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf(fp, "Level %d tmax %g. Ohf %3.2e, Ohd %3.2e, Ohs %3.2e, hf %3.2f\n", MAXlevel, tmax, Ohf, Ohd, Ohs, hf);
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
}
