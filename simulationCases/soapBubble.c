/**
 * @file wrinkling_bub_axi_v1.c
 * @author Saumili Jana (jsaumili@gmail.com)
 * @date 18-10-2024
 * Newtonian cases
 * 
 * Last update: Oct 28, 2024, Vatsal
 * changelog: fixed the initial condition. 
 * 
 * Last update: Mar 6, 2025, Saumili
 * Changelog: fixed boundary condition
 * 
 * Last update: Mar 15, 2025, Saumili
 * Changelog: updated pressure boundary conditions
 * 
 * Last update: Mar 31, 2025, Saumili
 *  Chhangelog: increased no of iterations for convergence
 * 
 *  * Last update: Apr 4, 2025, Saumili
 *  Chhangelog: corrected pressure initial condition 
*/

//f: 1 is liq, 0 is gas phase
#include "axi.h" //remove for 3d case
#include "navier-stokes/centered.h"
// #define FILTERED 1 // Smear density and viscosity jumps
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"//gravity

// error tolerances //for AMR
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)
#define MINlevel 2 

//time-intervals for saving
#define tsnap (0.01)  //snapshot saving interval

/*Id1 :indicates the liquid film  media formimg the bubble
Id2: indicates the surrounding gas/fluid(Newtonian)
*/
#define Rho21 (1e-3)
//#define Mu21 (1e-2)
//Calculations
#define Xcent (0.0)
#define Ycent (0.0)

#define R2circle(x,y) (sq(x - Xcent) + sq(y - Ycent))

//Boundary conditions
//velocity //x-axis axisymmetric
u.t[left] = dirichlet(0.0);//wall
u.n[left] = dirichlet(0.0);//wall
f[left] = neumann(0.0); // this sets the contact angle to 90 degrees.

u.t[right] = neumann(0.0);//outlow
u.n[right] = neumann(0.0);//outflow
//p[right] = neumann(0.0);//pressure outflow
p[right] = dirichlet(0.0);

u.t[top] = neumann(0.0);//outlow
u.n[top] = neumann(0.0);//outflow
//p[top] = neumann(0.0);//pressure outflow
p[top] = dirichlet(0.0);

//declarations
int MAXlevel;
double tmax, Oh1, Bo, Ldomain, k, h;

int main(int argc, char const *argv[]){
  //assignments
  //NITERMAX = 500; //increased no of iterations for convergence during initial timesteps for some cases
  MAXlevel = 8; //max possible grid res
  tmax = 1.0;
  Ldomain = 2.4;
  

  Bo = 1e-1; //gravity
  Oh1 = 1e-2;//liq film Oh

  k = 1e1; //curvature R/h

  fprintf(ferr, "Level %d, tmax %g, Bo %g, Oh1 %3.2e, Lo %g\n", MAXlevel, tmax, Bo, Oh1, Ldomain);

  L0=Ldomain;
  X0=-1.0; Y0=0.;
  init_grid (1 << (4));
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; rho2 = Rho21;
  f.sigma = 1; //coeff of surface tension
  mu1 = Oh1; mu2 = 1e-4;// mu2 = Mu21*Oh1;
  G.x = -Bo; //gravity
  run();
}


//Initial condition// 
event init(t = 0){
  if(!restore (file = "dump")){
    float y_p, x_p, x1, x2;
    h = 1/k;
    y_p = 0.1;
    x1 = sqrt(sq(1.0-h)-sq(y_p));
    x2 = sqrt(1-sq(y_p));
    x_p = (x1+x2)/2;


    refine((R2circle(x,y) < 1.05) && (R2circle(x,y) > sq(0.98-h)) && (level < MAXlevel));

    vertex scalar phi[];
    foreach_vertex() {
      if (y < y_p && x > 0.0) {
        // Lower part - half circle
        phi[] = (sq(h/2) - (sq(x - x_p) + sq(y - y_p)));
      } else {
        // Upper part - spherical rim
        double r = sqrt(sq(x) + sq(y));
        double shell = min(1. - r, (r - (1. - h)));
        phi[] = shell;
      }
    }
    fractions (phi, f);
    //pressure initialize

    foreach(){
      //p[] = (R2circle(x,y)<1)?4:0; //initialize pressure
      if ((R2circle(x,y)<sq(1.0-h)))
      {  p[] = 4;
      }
      else if ((R2circle(x,y)<=1)&&((R2circle(x,y)>=sq(1.0-h))))
      {
         p[] = 2;
      }
      else
      {
        p[] = 0;
      }
      u.x[] = 0;
      u.y[] = 0;
    }

  }
}

//AMR
scalar KAPPA[];
event adapt(i++){
  curvature(f, KAPPA);
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MAXlevel-4);
  //unrefine(x>150.0);

}

//Outpts
//static
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  p.nodump = false; // dump pressure
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {

  static FILE * fp;
  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t\n");
      fp = fopen ("log", "w");
      fprintf (fp, "Level %d, tmax %g, Oh %3.2e, Bo %2.1e, Lo %g\n", MAXlevel, tmax, Oh1, Bo, Ldomain);
      fprintf (fp, "i dt t\n");
      fprintf (fp, "%d %g %g \n", i, dt, t);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g\n", i, dt, t);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g\n", i, dt, t);
  }

}