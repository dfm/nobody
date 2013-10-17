#include "nobody.h"

#define INIT_ZERO(ntot,arr) { \
    for (int LOOP=0; LOOP < ntot; ++LOOP) arr[LOOP] = 0.0; \
}

#define DELTA_SUM_3(arr,i,j,delta,d2) { \
    int i3 = 3*i, j3 = j*3; \
    double value = arr[i3] - arr[j3]; \
    d2 = value * value; \
    delta[0] = value; \
    value = arr[i3+1] - arr[j3+1]; \
    d2 += value * value; \
    delta[1] = value; \
    value = arr[i3+2] - arr[j3+2]; \
    d2 += value * value; \
    delta[2] = value; \
}

#define SUM_SCALE_3(xin,xout,factor) { \
    double *XOUT = xout; \
    XOUT[0] += factor * xin[0]; \
    XOUT[1] += factor * xin[1]; \
    XOUT[2] += factor * xin[2]; \
}

nobody_system *nobody_allocate_system (int nbodies, double *masses)
{
    nobody_system *nbs = malloc(sizeof(nobody_system));
    nbs->nbodies = nbodies;
    nbs->masses = malloc(nbodies * sizeof(double));
    BLAS_dcopy(nbodies, masses, 1, nbs->masses, 1);
    return nbs;
}

void nobody_free_system (nobody_system *nbs)
{
    free(nbs->masses);
    free(nbs);
}

void nobody_gradient (double *x, double *aout, nobody_system *nbs)
{
    int i, j, n = nbs->nbodies, ntot = 3*n;
    double d2, delta[3], f1, f2, *m = nbs->masses;
    INIT_ZERO(ntot, aout);
    for (i = 0; i < n; ++i) {
        for (j = i+1; j < n; ++j) {
            DELTA_SUM_3 (x, i, j, delta, d2);
            d2 = d2 * sqrt(d2);
            f1 = -m[j] / d2;
            f2 = m[i] / d2;
            SUM_SCALE_3(delta, aout+3*i, f1);
            SUM_SCALE_3(delta, aout+3*j, f2);
        }
    }
}

void nobody_leapfrog_step (double h, double *x, double *v,
                           double *xout, double *vout, nobody_system *nbs)
{
    int i, N = 3*nbs->nbodies;
    double hh = 0.5 * h,
           *a = malloc(N*sizeof(double));

    // Update positions.
    nobody_gradient (x, a, nbs);
    for (i = 0; i < N; ++i) {
        vout[i] = v[i] + hh * a[i];
        xout[i] = x[i] + h * vout[i];
    }

    // Update velocities.
    nobody_gradient (xout, a, nbs);
    for (i = 0; i < N; ++i)
        vout[i] += hh * a[i];

    free(a);
}

void nobody_midpoint_step (double H, int nstep, double *x, double *v,
                        double *xout, double *vout, nobody_system *nbs)
{
    int i, k, n = 3*nbs->nbodies;
    double h = H / nstep, swap,
           *x0 = malloc(n * sizeof(double)),
           *x1 = malloc(n * sizeof(double)),
           *dx = malloc(n * sizeof(double)),
           *v0 = malloc(n * sizeof(double)),
           *v1 = malloc(n * sizeof(double)),
           *dv = malloc(n * sizeof(double));

    // Compute the initial gradient and save the first two terms in the
    // stencil using forward Euler.
    nobody_gradient (x, dv, nbs);
    for (i = 0; i < n; ++i) {
        x0[i] = x[i];
        x1[i] = x[i] + h * v[i];
        v0[i] = v[i];
        v1[i] = v[i] + h * dv[i];
    }

    // Loop over the requested number of midpoint steps and propagate the
    // gradient using centered difference.
    for (k = 2; k <= nstep; ++k) {
        nobody_gradient (x1, dv, nbs);
        for (i = 0; i < n; ++i) {
            swap = x1[i];
            x1[i] = x0[i] + 2 * h * v1[i];
            x0[i] = swap;

            swap = v1[i];
            v1[i] = v0[i] + 2 * h * dv[i];
            v0[i] = swap;
        }
    }

    // Compute the final position based on the last centered difference.
    nobody_gradient (x1, dv, nbs);
    for (i = 0; i < n; ++i) {
        xout[i] = 0.5 * (x0[i] + x1[i] + h * v1[i]);
        vout[i] = 0.5 * (v0[i] + v1[i] + h * dv[i]);
    }

    // Clean up.
    free(x0);
    free(x1);
    free(dx);
    free(v0);
    free(v1);
    free(dv);
}
