#ifndef _NOBODY_H_
#define _NOBODY_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.4625385377644

//
// Tools for solving Kepler's equation.
//
double kepler_ecc_to_mean_anomaly (double psi, double e);
double kepler_mean_to_ecc_anomaly (double M, double e, int *info);
int kepler_solve_one (int n, double *t, double mstar, double mplanet,
                      double e, double a, double t0, double pomega, double ix,
                      double iy, double *coords);
int kepler_solve (int n, double *t, double mstar, int np, double *m,
                  double *e, double *a, double *t0, double *pomega,
                  double *ix, double *iy, double *coords);

//
// Tools for solving a full N-body gravitational system.
//
void nobody_setup (int n, double *m);
void nobody_cleanup ();
void nobody_gradient (double t, double *x, double *xout);
double nobody_bs (int n, double t, double *x, double H, double eps,
                  double *xfinal, void (*dxdt) (double, double*, double*));

#endif
// _NOBODY_H_
