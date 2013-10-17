#include "nobody.h"

int main ()
{
    int nt = 10000, i, k, nbodies = 2;
    double masses[] = {1.0, 0.01},
           a[] = {100.0},
           t0[] = {0.0},
           e[] = {0.0},
           pomega[] = {0.0},
           ix[] = {0.0},
           iy[] = {0.0},
           x[6],
           v[6],
           *t = malloc(nt*sizeof(double)),
           *x_all = malloc(6*nt*sizeof(double)),
           *v_all = malloc(6*nt*sizeof(double)),
           h = 1e-4;

    for (i = 0; i < nt; ++i) t[i] = i * h;

    nobody_system *nbs = nobody_allocate_system (nbodies, masses);
    kepler_solve (nt, t, masses[0], 1, masses+1, a, t0, e, pomega, ix, iy,
                  x_all, v_all);

    for (i = 0; i < 6; ++i) {
        x[i] = x_all[i];
        v[i] = v_all[i];
    }

    for (i = 0; i < nt; ++i) {
        printf("%f %f %f %f %f %f\n", x[3]-x[0], x[4]-x[1], x[5]-x[2], x_all[i*6+3], x_all[i*6+4], x_all[i*6+5]);
        nobody_midpoint_step (h, 10, x, v, x, v, nbs);
    }

    nobody_free_system (nbs);
    free(t);
    free(x_all);
    free(v_all);

    return 0;
}
