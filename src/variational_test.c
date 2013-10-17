#include "nbody.h"

#define EQ_TOL        1e-4
#define N_TIME_KEPLER 5000

int main ()
{
    int i, j, k, info, nkick = 2;
    double t[N_TIME_KEPLER],
           coords_base[12*N_TIME_KEPLER],
           coords_var[12*N_TIME_KEPLER],
           mstar = 1.0,
           mplanet[1] = {1.0e-5},
           e[1]       = {0.4},
           a[1]       = {150.0},
           t0[1]      = {10.0},
           pomega[1]  = {-0.5},
           ix[1]      = {0.5},
           iy[1]      = {1.6},
           kicks[2]   = {25.0, 100.0},
           del_e[3]   = {0.0, 0.0, 0.0},
           del_a[3]   = {0.0, 0.0, 0.0},
           del_t0[3]  = {0.0, 0.0, 0.0},
           del_pomega[3] = {0.0, 0.0, 0.0},
           del_ix[3]  = {0.0, 0.0, 0.0},
           del_iy[3]  = {0.0, 0.0, 0.0},
           delta;

    // Initialize the time array.
    for (k = 0; k < N_TIME_KEPLER; ++k)
        t[k] = k * 200.0 / N_TIME_KEPLER;

    // Solve the pure Kepler system.
    info = kepler_solve (N_TIME_KEPLER, t, mstar, 1, mplanet, a, t0, e,
                         pomega, ix, iy, coords_base);
    if (info != 0) {
        fprintf(stderr, "Base Kepler solve failed in variational test.\n");
        return 1;
    }

    info = variational_solve (N_TIME_KEPLER, t, mstar, 1, mplanet, a, t0,
                              e, pomega, ix, iy,
                              nkick, kicks, del_a, del_t0, del_e,
                              del_pomega, del_ix, del_iy, coords_var);
    if (info != 0) {
        fprintf(stderr, "Variational Kepler solve failed.\n");
        return 1;
    }

    // Check that all of the values are equal.
    for (k = 0; k < N_TIME_KEPLER; ++k) {
        for (i = 0; i < 12; ++i) {
            delta = coords_base[12*k+i] - coords_var[12*k+i];
            delta = sqrt(delta * delta);
            if (delta > EQ_TOL) {
                fprintf(stderr,
                        "Variational and base coordinates are not equal.\n");
                fprintf(stderr, "%d %d %e %e %e\n",
                        k, i, t[k], coords_base[12*k+i], coords_var[12*k+i]);
                return 1;
            }
        }
    }
    return 0;
}
