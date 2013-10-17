#include "nbody.h"

#define EQ_TOL        1e-4

int main (int argc, char *argv[])
{
    int info;
    double M, psi_est, delta, psi, e, fail;

    if (argc < 3) {
        psi = -10.0;
        e = 0.9;
        fail = 0;
    } else {
        psi = atof(argv[1]);
        e = atof(argv[2]);
        if (argc < 4) fail = 0;
        else fail = atoi(argv[3]);
    }

    M = kepler_ecc_to_mean_anomaly(psi, e);
    psi_est = kepler_mean_to_ecc_anomaly(M, e, &info);
    if (info != 0) {
        if (info == fail) return 0;
        fprintf(stderr, "Solver failed for e=%e and psi=%e\n", e, psi);
        return 1;
    }

    delta = fmod(psi, 2 * M_PI) - psi_est;
    delta = sqrt(delta * delta);
    if (delta <= EQ_TOL) return 0;

    fprintf(stderr,
            "Solver not close enough for e=%e and psi=%e. psi_est=%e\n",
            e, psi, psi_est);
    return 1;
}
