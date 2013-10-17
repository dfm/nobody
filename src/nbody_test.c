#include "nbody.h"

#define EQ_TOL        1e-4

int main ()
{
    int i, j, k, info;
    double t[1] = {0.0},
           coords[18],
           tmp[12],
           mstar = 1.0,
           mplanet[1] = {1.0e-5},
           e[1]       = {0.4},
           a[1]       = {150.0},
           t0[1]      = {0.0},
           pomega[1]  = {0.0},
           ix[1]      = {0.0},
           iy[1]      = {0.0};

    int N = 5000;
    double initial[12], masses[2] = {mstar, mplanet[0]},
           time = t[0], dt = 0.5, diff, maxdiff;

    // Set up N-body initial conditions.
    nbody_system *nbody = nbody_setup(2, masses);
    info = kepler_solve_one (1, t, mstar, mplanet[0], a[0], t0[0], e[0],
                             pomega[0], ix[0], iy[0], coords);

    for (i = 0; i < 6; ++i) {
        initial[i] = 0.0;
        initial[6+i] = coords[i];
    }

    // Loop over time steps.
    for (k = 0; k < N; ++k) {
        // Advance the time using the N-body simulation.
        time = nbody_bs (12, time, initial, dt, 1e-12, initial,
                         nbody, &info);

        // Solve the Kepler version of the problem.
        t[0] = time;
        kepler_solve (1, t, mstar, 1, mplanet, a, t0, e, pomega, ix, iy, tmp);

        // Compute the differences.
        maxdiff = 0.0;
        for (i = 0; i < 6; ++i) {
            diff = tmp[6+i] - (initial[6+i] - initial[i]);
            diff = sqrt(diff * diff);
            if (diff > maxdiff) maxdiff = diff;
        }

        // Check for equality.
        if (maxdiff >= EQ_TOL) {
            fprintf(stderr, "Orbit position off by too much at t=%e\n",
                    time);
            fprintf(stderr, "%e\n", maxdiff);
            return 1;
        }
    }

    return 0;
}
