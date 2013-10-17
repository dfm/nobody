#include "nbody.h"

#define EQ_TOL        1e-4

void test_one_velocity (int *total, int *passed, double time,
                        double e, double pomega, double t0, double ix,
                        double iy)
{
    int i, j, info;
    double eps = 1e-4,
           t[3] = {time-eps, time, time+eps},
           coords[18],
           mstar = 1.0,
           mplanet = 1.0e-5,
           a = 150.0,
           numerical, analytic, delta;

    info = kepler_solve_one (3, t, mstar, mplanet, a, t0, e, pomega,
                             ix, iy, coords);
    if (info != 0) return;

    for (i = 0; i < 3; ++i) {
        numerical = 0.5 * (coords[12+i] - coords[i]) / eps;
        analytic  = coords[9+i];
        delta = numerical - analytic;
        delta *= delta;
        delta = sqrt(delta);

        if (delta <= EQ_TOL) ++(*passed);
        else {fprintf(stderr,
                     "Error in velocity %d too large (%e).\n",
                     i, delta);
            printf("pomega = %e, %e %e\n", pomega, numerical, analytic);
        }

        ++(*total);
    }
}

int main ()
{
    int total = 0, passed = 0;
    double e, pomega, t, t0, ix, iy;
    for (e = 0.0; e < 1.0; e += 0.1)
    for (pomega = -6.0; pomega < 6.0; pomega += 2.0)
    for (t0 = -10.0; t0 < 10.0; t0 += 2.0)
    for (t = 0.0; t < 1000.0; t += 50.0)
    for (ix = -6.0; ix < 6.0; ix += 2.0)
    for (iy = -6.0; iy < 6.0; iy += 2.0)
        test_one_velocity(&total, &passed, t, e, pomega, t0, ix, iy);
    if (total != passed) return 1;
    return 0;
}
