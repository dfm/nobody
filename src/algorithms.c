#include <math.h>
#include <stdlib.h>

double compute_acceleration (int nbody, double *x, double *m)
{
    int i, j, k;

    double dr, r, r2, r3, dv, v2, da, a2, dm, epot = 0.0;

    // Update the acceleration.
    for (i = 0; i < nbody; ++i)
        for (k = 0; k < 3; ++k)
            x[9 * i + k + 6] = 0.0;
    for (i = 0; i < nbody; ++i) {
        for (j = i + 1; j < nbody; ++j) {
            // Compute the phase space distance between particles.
            r2 = 0.0;
            v2 = 0.0;
            for (k = 0; k < 3; ++k) {
                dr = x[9 * i + k] - x[9 * j + k];
                dv = x[9 * i + k + 3] - x[9 * j + k + 3];
                r2 += dr * dr;
                v2 += dv * dv;
            }
            r = sqrt(r2);
            r3 = r * r2;

            // Update the potential energy.
            epot -= m[i] * m[j] / r;

            // Compute the acceleration.
            a2 = 0.0;
            for (k = 0; k < 3; ++k) {
                da = (x[9 * j + k] - x[9 * i + k]) / r3;
                x[9 * i + k + 6] += da * m[j];
                x[9 * j + k + 6] -= da * m[i];
            }
        }
    }

    return epot;
}

double euler (int nbody, double *x, double *m, double dt)
{
    int i, k;
    double ekin = 0.0, epot;

    epot = compute_acceleration(nbody, x, m);

    // Update the position and velocity.
    for (i = 0; i < nbody; ++i) {
        for (k = 0; k < 6; ++k)
            x[9 * i + k] += x[9 * i + k + 3] * dt;

        // Update the kinetic energy.
        for (k = 0; k < 3; ++k)
            ekin += m[i] * x[9 * i + 3 + k] * x[9 * i + 3 + k];
    }

    return ekin + epot;
}

double leapfrog (int nbody, double *x, double *m, double dt)
{
    int i, k;
    double ekin = 0.0, epot;

    // Initial half step.
    for (i = 0; i < nbody; ++i)
        for (k = 0; k < 3; ++k) {
            x[9 * i + k + 3] += 0.5 * x[9 * i + k + 6] * dt;
            x[9 * i + k] += x[9 * i + k + 3] * dt;
        }

    epot = compute_acceleration(nbody, x, m);

    // Update the velocity.
    for (i = 0; i < nbody; ++i)
        for (k = 3; k < 6; ++k) {
            x[9 * i + k] += 0.5 * x[9 * i + k + 3] * dt;
            ekin += m[i] * x[9 * i + k] * x[9 * i + k];
        }

    return ekin + epot;
}
