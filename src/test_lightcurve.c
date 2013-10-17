#import "nbody.h"
#import "lightcurve.h"

#define N_TIME_KEPLER 5000
#define N_PLANETS     1
#define N_LDP         100

int main ()
{
    int i, j, k, info, nkick = 2;
    double t[N_TIME_KEPLER],
           coords[12*N_TIME_KEPLER],
           mstar = 1.0,
           rstar = 1.0,
           mplanet[N_PLANETS] = {1.0e-5},
           rp[N_PLANETS]      = {0.1},
           e[N_PLANETS]       = {0.9},
           a[N_PLANETS]       = {10.0},
           t0[N_PLANETS]      = {1.0},
           pomega[N_PLANETS]  = {1.5 * M_PI},
           ix[N_PLANETS]      = {0.0},
           iy[N_PLANETS]      = {0.0},
           rld[N_LDP],
           ild[N_LDP],
           flux[N_TIME_KEPLER],
           fp[N_TIME_KEPLER];

    // Set up the limb darkening.
    for (i = 0; i < N_LDP; ++i) {
        rld[i] = (i+1) * 1.0 / N_LDP;
        double omm = 1 - sqrt(1-rld[i]*rld[i]);
        ild[i] = 1.0 - 0.3 * omm - 0.1 * omm * omm;
    }

    // Initialize the time array.
    for (k = 0; k < N_TIME_KEPLER; ++k) {
        t[k] = k * 10.0 / N_TIME_KEPLER;
        fp[k] = 1e-4;
    }

    // Solve the Kepler system.
    info = kepler_solve (N_TIME_KEPLER, t, mstar, N_PLANETS, mplanet, a,
                         t0, e, pomega, ix, iy, coords);

    // Compute the light curve.
    lightcurve_full (N_TIME_KEPLER, N_PLANETS, coords, fp, rstar, rp, N_LDP,
                     rld, ild, flux);

    /* for (i = 0; i < N_TIME_KEPLER; ++i) */
    /*     printf("%e %e %e %e\n", t[i], flux[i], coords[12*i+6], coords[12*i+7]); */

    return 0;
}
