#include "nobody.h"

#define EQ_TOL        1e-4
#define N_TIME_KEPLER 5000

int test_one_solver (double psi, double e, int fail)
{
    int info;
    double M, psi_est, delta;

    M = kepler_ecc_to_mean_anomaly(psi, e);
    psi_est = kepler_mean_to_ecc_anomaly(M, e, &info);
    if (info != 0) {
        if (info == fail) return 1;
        fprintf(stderr, "Solver failed for e=%e and psi=%e\n", e, psi);
        return 0;
    }

    delta = fmod(psi, 2 * M_PI) - psi_est;
    delta = sqrt(delta * delta);
    if (delta <= EQ_TOL) return 1;

    fprintf(stderr,
            "Solver not close enough for e=%e and psi=%e. psi_est=%e\n",
            e, psi, psi_est);
    return 0;
}

void test_solver (int *total, int *passed)
{
    double e;

    *passed += test_one_solver(-10.0, 0.9, 0);
    ++(*total);

    *passed += test_one_solver(100.0, 0.8, 0);
    ++(*total);

    *passed += test_one_solver(100.0, 1.0, 2);
    ++(*total);

    for (e = 1.0 - 1e-10; e < 1.0; e += 1e-11) {
        *passed += test_one_solver(1.0, e, 0);
        ++(*total);
    }

    *passed += test_one_solver(M_PI, 0.6, 0);
    ++(*total);

    *passed += test_one_solver(0.5 * M_PI, 0.6, 0);
    ++(*total);

    *passed += test_one_solver(2 * M_PI, 0.9, 0);
    ++(*total);
}

void test_velocities (int *total, int *passed)
{
    int i, j, info;
    double eps = 1e-4,
           t[3] = {0.0, eps, 2 * eps},
           coords[18],
           mstar = 1.0,
           mplanet = 1.0e-5,
           e = 0.0, a = 150.0, t0 = 0.0, pomega = 0.5, ix = 0.1, iy = 0.1,
           numerical, analytic, delta;

    info = kepler_solve_one (3, t, coords, mstar, mplanet, e, a, t0, pomega,
                             ix, iy);

    for (i = 0; i < 3; ++i) {
        numerical = 0.5 * (coords[6*2+i] - coords[i]) / eps;
        analytic  = coords[9+i];
        delta = (numerical - analytic) / numerical;
        delta *= delta;
        delta = sqrt(delta);

        if (delta <= EQ_TOL) ++(*passed);
        else fprintf(stderr,
                     "Relative error in velocity %d is not acceptable (%e).\n",
                     i, delta);

        ++(*total);
    }
}

void test_kepler (int *total, int *passed)
{
    int i, j, info;
    double t[N_TIME_KEPLER],
           coords[6 * N_TIME_KEPLER],
           mstar = 1.0,
           mplanet = 1.0e-5,
           e = 0.4, a = 150.0, t0 = 0.0, pomega = 0.0, ix = 0.0, iy = 0.0;

    // Set up the time array.
    for (i = 0; i < N_TIME_KEPLER; ++i)
        t[i] = 1000. * i / N_TIME_KEPLER;

    // Solve the system.
    info = kepler_solve_one (N_TIME_KEPLER, t, coords, mstar, mplanet, e, a,
                             t0, pomega, ix, iy);
}

void test_nbody (int *total, int *passed)
{
    int i, j, k, info;
    double eps = 1e-4,
           t[3] = {0.0, eps, 2 * eps},
           coords[18],
           mstar = 1.0,
           mplanet = 1.0e-5,
           e = 0.4, a = 150.0, t0 = 0.0, pomega = 0.0, ix = 0.0, iy = 0.0;

    int N = 1000;
    double initial[12], masses[2] = {mstar, mplanet}, time = t[1], dt = 0.5;

    // Set up N-body initial conditions.
    nobody_setup(2, masses);
    info = kepler_solve_one (3, t, coords, mstar, mplanet, e, a, t0, pomega,
                             ix, iy);

    for (i = 0; i < 3; ++i) {
        // Planet position.
        initial[6+i] = coords[6+i];

        // Planet velocity.
        initial[9+i] = 0.5 * (coords[6*2+i] - coords[i]) / eps;

        // Star coordinates.
        initial[i] = 0.0;
        initial[3+i] = 0.0;
    }

    // Loop over time steps.
    for (k = 0; k < N; ++k) {
        t[0] = time;
        info = kepler_solve_one (1, t, coords, mstar, mplanet, e, a, t0,
                                 pomega, ix, iy);

        time = nobody_bs (12, time, initial, dt, 1e-9, initial,
                          &nobody_gradient);

        printf("%e ", time);
        for (i = 0; i < 12; ++i) printf("%e ", initial[i]);
        for (i = 0; i < 6; ++i) printf("%e ", coords[i]);
        printf("\n");
    }
}

int main ()
{
    int total = 0, passed = 0;

    test_solver (&total, &passed);
    test_velocities (&total, &passed);
    test_kepler (&total, &passed);
    test_nbody (&total, &passed);

    fprintf(stderr, "%d test(s) passed out of %d tried.\n", passed, total);

    return 0;
}
