#ifndef _NOBODY_H_
#define _NOBODY_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nobody_blas.h"

//
// Newton's constant in $R_\odot^3 M_\odot^{-1} {days}^{-2}$.
//
#define G_GRAV 2945.4625385377644

//
// Tools for solving Kepler's equation.
//
double kepler_ecc_to_mean_anomaly (double psi, double e);
double kepler_mean_to_ecc_anomaly (double M, double e, int *info);
int kepler_solve_one (
    int n,
    double *t,
    double mstar,
    double mplanet,
    double a,
    double t0,
    double e,
    double pomega,
    double ix,
    double iy,
    double *x,
    double *v,
    int incr
);
int kepler_solve (
    int n,
    double *t,
    double mstar,
    int np,
    double *m,
    double *a,
    double *t0,
    double *e,
    double *pomega,
    double *ix,
    double *iy,
    double *x,
    double *v
);

//
// Tools for solving a variational Kepler system.
//
int variational_solve (int n, double *t, double mstar, int np,
                       double *m, double *a, double *t0, double *e,
                       double *pomega, double *ix, double *iy,
                       int nkick, double *kicks,
                       double *del_a, double *del_t0, double *del_e,
                       double *del_pomega, double *del_ix, double *del_iy,
                       double *coords);

//
// Tools for solving a full N-body gravitational system.
//

typedef struct nobody_system_struct {

    int nbodies;
    double *masses;

} nobody_system;

nobody_system *nobody_allocate_system (
    int nbodies,                // The number of bodies in the system.
    double *masses              // The masses of the bodies (len==nbodies).
);

void nobody_free_system (
    nobody_system *nbs          // The system to free.
);

void nobody_gradient (
    double *x,                  // The current positions of the bodies.
    double *aout,               // The instantaneous accelerations.
    nobody_system *nbs          // The current nobody workspace.
);

void nobody_leapfrog_step (
    double h,
    double *x,
    double *v,
    double *xout,
    double *vout,
    nobody_system *nbs
);
void nobody_midpoint_step (double H, int nstep, double *x, double *v,
                        double *xout, double *vout, nobody_system *nbs);

#endif
// _NOBODY_H_
