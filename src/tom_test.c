#include "nbody.h"
#include <stdio.h>

int main ()
{
    int i, j, info;
    double mp[2] = {0.01, 0.01},
           masses[3] = {1.0, mp[0], mp[1]},
           a[2] = {4.2238165436989386, 6.926564935838537},
           t0[2] = {0.0, 0.0},
           zeros[2] = {0.0, 0.0},
           t, time = 0.0, ta[1] = {0.0},
           initial[18];

    nbody_setup(3, masses);

    // Set up the initial conditions.
    info = kepler_solve (1, ta, masses[0], 2, mp, a, t0, zeros, zeros,
                         zeros, zeros, initial);
    if (info) {
        printf("Failed to initialize.\n");
        return info;
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 6; ++j)
            printf("%e\t", initial[i*6+j]);
        printf("\n");
    }

    for (t = 0.02; t < 30.0; t += 0.02) {
        time = nbody_bs (18, time, initial, t - time, 1e-12, initial,
                         &nbody_gradient, &info);
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 6; ++j)
                printf("%e\t", initial[i*6+j]);
            printf("\n");
        }
    }

    printf("%f\n", time);

    return 0;
}
