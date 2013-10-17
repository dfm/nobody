#include "nbody.h"

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

void test_velocities (int *total, int *passed)
{
    double e, pomega, t, t0, ix, iy;
    for (e = 0.0; e < 1.0; e += 0.1)
    for (pomega = -6.0; pomega < 6.0; pomega += 2.0)
    for (t0 = -10.0; t0 < 10.0; t0 += 2.0)
    for (t = 0.0; t < 1000.0; t += 50.0)
    for (ix = -6.0; ix < 6.0; ix += 2.0)
    for (iy = -6.0; iy < 6.0; iy += 2.0)
        test_one_velocity(total, passed, t, e, pomega, t0, ix, iy);
}

void test_nbody (int *total, int *passed)
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

    int N = 1000;
    double initial[12], masses[2] = {mstar, mplanet[0]},
           time = t[0], dt = 0.5, diff, maxdiff;

    ++(*total);

    // Set up N-body initial conditions.
    nbody_setup(2, masses);
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
                          &nbody_gradient, &info);

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
            return;
        }
    }
    ++(*passed);
}

void test_variational (int *total, int *passed)
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

    ++(*total);

    // Initialize the time array.
    for (k = 0; k < N_TIME_KEPLER; ++k)
        t[k] = k * 200.0 / N_TIME_KEPLER;

    // Solve the pure Kepler system.
    info = kepler_solve (N_TIME_KEPLER, t, mstar, 1, mplanet, a, t0, e,
                         pomega, ix, iy, coords_base);
    if (info != 0) {
        fprintf(stderr, "Base Kepler solve failed in variational test.\n");
        return;
    }

    info = variational_solve (N_TIME_KEPLER, t, mstar, 1, mplanet, a, t0,
                              e, pomega, ix, iy,
                              nkick, kicks, del_a, del_t0, del_e,
                              del_pomega, del_ix, del_iy, coords_var);
    if (info != 0) {
        fprintf(stderr, "Variational Kepler solve failed.\n");
        return;
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
                return;
            }
        }
    }

    ++(*passed);
}

int main ()
{
    int total = 0, passed = 0;

    test_solver (&total, &passed);
    test_velocities (&total, &passed);
    test_nbody (&total, &passed);
    test_variational (&total, &passed);

    fprintf(stderr, "%d test(s) passed out of %d tried.\n", passed, total);

    return 0;
}
