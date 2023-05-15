/**
# Three-phase interfacial flows

This file helps setup simulations for flows of three fluids separated by
corresponding interfaces (i.e. immiscible fluids). It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction in fluid i is $f_i=1$ and $f_i=0$ everywhere else.
The densities and dynamic viscosities for fluid 1, 2 and 3 are *rho1*,
*mu1*, *rho2*, *mu2*, and *rho3*, *mu3* respectively. */

#include "vof.h"
scalar f1[], f2[], *interfaces = {f1, f2};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;

/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */

  if (mu1 || mu2)
    mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
#define rho(f1, f2) (clamp(f1*(1-f2), 0., 1.) * rho1 + clamp(f1*f2, 0., 1.) * rho2 + clamp((1-f1), 0., 1.) * rho3)
#endif
#ifndef mu
#define mu(f1, f2) (clamp(f1*(1-f2), 0., 1.) * mu1 + clamp(f1*f2, 0., 1.) * mu2 + clamp((1-f1), 0., 1.) * mu3)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event tracer_advection (i++)
{
  if (i > 1){
  foreach(){
    if ((f2[] > 0.5) && (f1[] < 0.5)){
      f1[] = f2[];
    }
  }
}

/**
When using smearing of the density jump, we initialise *sfi* with the
vertex-average of *fi*. */
#ifdef FILTERED
  int counter1 = 0;
  for (scalar sf in smearInterfaces){
    counter1++;
    int counter2 = 0;
    for (scalar f in interfaces){
      counter2++;
      if (counter1 == counter2){
        // fprintf(ferr, "%s %s\n", sf.name, f.name);
      #if dimension <= 2
          foreach(){
            sf[] = (4.*f[] +
        	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
        	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
          }
      #else // dimension == 3
          foreach(){
            sf[] = (8.*f[] +
        	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
        	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
        		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
        		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
        	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
        	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
          }
      #endif
      }
    }
  }
#endif

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = refine_bilinear;
    sf.dirty = true; // boundary conditions need to be updated
  }
#endif
}


event properties (i++) {
  
  foreach_face() {
  double ff1 = (sf1[] + sf1[-1])/2.;
  double ff2 = (sf2[] + sf2[-1])/2.;
  /**
  This is a bit of fine tuning here. We might get places in the domain
  where density and viscosity go below those of air, we just make the properties
  at those locations equal to that of air
  */
  alphav.x[] = fm.x[]/rho(ff1, ff2);
  face vector muv = mu;
  muv.x[] = fm.x[]*mu(ff1, ff2);
  }

  foreach(){
  rhov[] = cm[]*rho(sf1[], sf2[]);
  }

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    sf.dirty = true;
  }
#endif
}
