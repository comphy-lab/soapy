/**
# Axisymmetric Bubble Wrinkling Simulation

## Overview
This simulation models the dynamics of a thin liquid film forming a bubble
with wrinkling instabilities. The code solves the two-phase Navier-Stokes
equations in an axisymmetric configuration, capturing the interface dynamics
between a liquid film and surrounding gas phase.

## Physics Description
The simulation captures:
- Two-phase flow dynamics with surface tension effects
- Gravitational forcing characterized by Bond number (Bo)
- Viscous effects characterized by Ohnesorge number (Oh)
- Interface tracking using Volume-of-Fluid (VOF) method
- Adaptive mesh refinement for efficient computation

## Authors and Changelog
- **Author**: Saumili Jana (jsaumili@gmail.com)
- **Created**: October 18, 2024

### Version History
- **Oct 28, 2024** (Vatsal): Fixed initial condition implementation
- **Mar 6, 2025** (Saumili): Corrected boundary condition formulation
- **Mar 15, 2025** (Saumili): Updated pressure boundary conditions
- **Mar 31, 2025** (Saumili): Increased iterations for convergence
- **Apr 4, 2025** (Saumili): Corrected pressure initial condition

## Phase Convention
- f = 1: Liquid phase (thin film)
- f = 0: Gas phase (surrounding medium)
*/

#include "axi.h"  // Axisymmetric coordinate system
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
scalar T[];
scalar * tracers = {T};
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // For gravitational acceleration

/**
## Simulation Parameters

### Adaptive Mesh Refinement Tolerances
These control the error thresholds for grid adaptation:
- Volume fraction error tolerance
- Velocity component error tolerance
- Curvature error tolerance
- Acceleration error tolerance
*/
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)
#define MINlevel 2

/**
### Time Control Parameters
- `tsnap`: Interval for saving simulation snapshots
*/
#define tsnap (0.01)

/**
### Physical Properties
- `Rho21`: Density ratio between gas and liquid phases (ρ₂/ρ₁)
- `Xcent`, `Ycent`: Center coordinates of the bubble
- `R2circle`: Macro for computing squared radial distance
*/
#define Rho21 (1e-3)
#define Xcent (0.0)
#define Ycent (0.0)
#define R2circle(x,y) (sq(x - Xcent) + sq(y - Ycent))

/**
## Boundary Conditions

### Left Boundary (Axis of Symmetry)
- No-slip wall condition for velocity
- Contact angle set to 90 degrees via Dirichlet condition on f
*/
u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

/**
### Right Boundary (Outflow)
- Neumann conditions for velocity components
- Dirichlet condition for pressure (atmospheric)
*/
u.t[right] = neumann(0.0);
u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);

/**
### Top Boundary (Outflow)
- Neumann conditions for velocity components
- Dirichlet condition for pressure (atmospheric)
*/
u.t[top] = neumann(0.0);
u.n[top] = neumann(0.0);
p[top] = dirichlet(0.0);

/**
## Global Variables
- `MAXlevel`: Maximum refinement level for adaptive mesh
- `tmax`: Maximum simulation time
- `Oh1`: Ohnesorge number for liquid film
- `Bo`: Bond number (gravity/surface tension ratio)
- `Ldomain`: Domain size
- `k`: Curvature parameter (R/h ratio)
- `h`: Film thickness
*/
int MAXlevel;
double tmax, Oh1, Bo, Ldomain, k, h;

/**
## Main Function

### Purpose
Initializes simulation parameters and launches the computation.

### Key Parameters Set
- Maximum grid refinement level: 8
- Simulation duration: 1.0 time units
- Domain size: 2.4 units
- Bond number: 0.1 (weak gravity effects)
- Ohnesorge number: 0.01 (moderate viscous effects)
- Curvature parameter: 10 (thin film approximation)
*/
int main(int argc, char const *argv[]) {
  // Parameter assignments
  MAXlevel = 9;
  tmax = 1.0;
  Ldomain = 2.4;
  
  Bo = 1e-1;
  Oh1 = 1e-2;
  k = 1e1;

  fprintf(ferr, "Level %d, tmax %g, Bo %g, Oh1 %3.2e, Lo %g\n", 
          MAXlevel, tmax, Bo, Oh1, Ldomain);

  // Domain configuration
  L0 = Ldomain;
  X0 = -1.01;
  Y0 = 0.;
  init_grid(1 << 4);
  
  // Create output directory
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // Physical properties assignment
  rho1 = 1.0;
  rho2 = Rho21;
  f.sigma = 1;  // Surface tension coefficient
  mu1 = Oh1;
  mu2 = 1e-4;
  G.x = -Bo;    // Gravitational acceleration

  run();
}

/**
## Initial Condition Event

### Purpose
Sets up the initial bubble geometry with a half-circle lower part and
spherical rim upper part. This creates a bubble with controlled wrinkling
potential based on the curvature parameter k.

### Geometry Construction
- Lower region: Half-circular cap
- Upper region: Spherical shell rim
- Film thickness: h = 1/k
- Initial refinement focused near the interface
*/
event init(t = 0) {
  if (!restore(file = "dump")) {
    float y_p, x_p, x1, x2;
    h = 1/k;
    y_p = 0.1;
    x1 = sqrt(sq(1.0 - h) - sq(y_p));
    x2 = sqrt(1 - sq(y_p));
    x_p = (x1 + x2) / 2;

    // Adaptive refinement near interface
    refine((R2circle(x, y) < 1.05) && 
           (R2circle(x, y) > sq(0.98 - h)) && 
           (level < MAXlevel));

    // Initialize level-set function for interface
    vertex scalar phi[];
    foreach_vertex() {
      if (y < y_p && x > 0.0) {
        // Lower part - half circle geometry
        phi[] = (sq(h/2) - (sq(x - x_p) + sq(y - y_p)));
      } else {
        // Upper part - spherical rim
        double r = sqrt(sq(x) + sq(y));
        double shell = min(1. - r, (r - (1. - h)));
        phi[] = shell;
      }
    }
    fractions(phi, f);

    // Initialize pressure field based on region
    foreach() {
      if (R2circle(x, y) < sq(1.0 - h)) {
        p[] = 4;  // Inner bubble pressure
      }
      else if ((R2circle(x, y) <= 1) && 
               (R2circle(x, y) >= sq(1.0 - h))) {
        p[] = 2*f[];  // Film pressure
      }
      else {
        p[] = 0;  // Atmospheric pressure
      }
      u.x[] = 0;
      u.y[] = 0;
    }
    fraction(T, (sq(1.0 - h) - R2circle(x, y)));
  }
}

/**
## Adaptive Mesh Refinement Event

### Purpose
Dynamically refines the computational mesh based on solution gradients.
This ensures high resolution near the interface and regions of high
curvature while maintaining efficiency in smooth regions.

### Refinement Criteria
- Volume fraction gradients (interface tracking)
- Velocity component gradients (flow features)
- Curvature gradients (wrinkling detection)
*/
scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA},
                (double[]){fErr, VelErr, VelErr, KErr},
                MAXlevel, MAXlevel - 4);
}

/**
## Data Output Event

### Purpose
Saves simulation state at regular intervals for post-processing and
restart capability.

### Output Schedule
- Initial condition saved at t = 0
- Snapshots saved every tsnap time units
- Full field data including pressure dumped
*/
event writingFiles(t = 0, t += tsnap; t <= tmax) {
  p.nodump = false;  // Enable pressure output
  dump(file = "dump");
  
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}

/**
## Logging Event

### Purpose
Maintains a detailed log of simulation progress including timestep
information and key parameters.

### Log Contents
- Iteration number
- Current timestep size
- Simulation time
- Initial header with all simulation parameters
*/
event logWriting(i++) {
  static FILE * fp;
  if (pid() == 0) {
    if (i == 0) {
      fprintf(ferr, "i dt t\n");
      fp = fopen("log", "w");
      fprintf(fp, "Level %d, tmax %g, Oh %3.2e, Bo %2.1e, Lo %g\n", 
              MAXlevel, tmax, Oh1, Bo, Ldomain);
      fprintf(fp, "i dt t\n");
      fprintf(fp, "%d %g %g \n", i, dt, t);
      fclose(fp);
    } else {
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g\n", i, dt, t);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g\n", i, dt, t);
  }
}