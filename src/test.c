#include <stdio.h>
#include <stdlib.h>

double compute_acceleration (int nbody, double *x, double *m);
double euler (int nbody, double *x, double *m, double dt);
double leapfrog (int nbody, double *x, double *m, double dt);

int main ()
{
    int nbody = 2, i, k, n;
    double energy, dt = 1.256e-4;
    double *m = malloc(nbody * sizeof(double)),
           *x = malloc(9 * nbody * sizeof(double));

    for (i = 0; i < nbody * 9; ++i) x[i] = 0.0;
    x[0] = 1.0;
    x[4] = 0.5;
    // x[13] = -0.02;
    m[0] = 0.1;
    m[1] = 10.0;

    // Compute the initial acceleration.
    compute_acceleration(nbody, x, m);

    for (n = 0; n < 1e5; ++n) {
        energy = leapfrog(nbody, x, m, dt);
        for (i = 0; i < 9 * nbody; ++i) printf("%.10e ", x[i]);
        printf("%e\n", energy);
    }

    free(x);
    free(m);
    return 0;
}
