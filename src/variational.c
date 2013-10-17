#include "nbody.h"

//
// Solve Kepler's equation for the orbits of a system of planets.
//
// :param int n:             The number of time steps to solve.
// :param double t[n]:       The times where the orbit should be computed.
//                           This array is assumed to be sorted.
// :param double mstar:      The mass of the star in Solar masses.
// :param int np:            The number of planets in the system.
// :param double m[np]:      The masses of the planets in Solar masses.
// :param double a[np]:      The semi-major axis of the mean orbits in Solar
//                           radii.
// :param double t0[np]:     The epochs of the mean orbits (defined as the
//                           center time of a target transit) in days.
// :param double e[np]:      The eccentricities of the mean orbits.
// :param double pomega[np]: The orientations of the mean orbital ellipses.
//                           Defined as the angle between the major-axis of
//                           the ellipse and the observer in radians.
// :param double ix[np]:     The inclinations of the mean orbits relative to
//                           the observer in radians. At ``ix=0``, the system
//                           is edge on and at ``ix=0.5*M_PI``, the system is
//                           face on.
// :param double iy[np]:     The orientations of the mean orbits in the plane
//                           of the sky measured in radians.
// :param int nkick:         The number of kicks that the system will receive.
// :param double kicks[nkick]:
//                           The times of the kicks in days. This array is
//                           assumed to be sorted.
// :param double del_a[(nkick+1)*np]:
//                           The change in the semi-major axes of the orbits
//                           after each kick.
// :param double del_t0[(nkick+1)*np]:
//                           The change in the epochs of the orbits after each
//                           kick.
// :param double del_e[(nkick+1)*np]:
//                           The change in the eccentricities of the orbits
//                           after each kick.
// :param double del_pomega[(nkick+1)*np]:
//                           The change in the angle of the orbital ellipses
//                           after each kick.
// :param double del_ix[(nkick+1)*np]:
//                           The change in the inclinations of the orbital
//                           planes after each kick.
// :param double del_iy[(nkick+1)*np]:
//                           The change in the inclinations of the orbital
//                           planes in the plane of the sky after each kick.
// :param double coords[6*(np+1)*n]:
//                           The phase space coordinates of the star and the
//                           planets at each time ``t``. The coordinate system
//                           is such that the first element points towards the
//                           observer and the other two are in the plane of
//                           the sky (as defined by ``iy``). The positions
//                           are measured in Solar radii and the velocities
//                           in Solar radii per day.
//
// :returns int info:        A flag describing the success of the computation.
//                           0=success, 1=solve did not converge after 500
//                           iterations, 2=unphysical eccentricity.
//
int variational_solve (int n, double *t, double mstar, int np,
                       double *m, double *a, double *t0, double *e,
                       double *pomega, double *ix, double *iy,
                       int nkick, double *kicks,
                       double *del_a, double *del_t0, double *del_e,
                       double *del_pomega, double *del_ix, double *del_iy,
                       double *coords)
{
    int i, kick, ind = 0, start, count, exitflag, info = 0;
    double *ttmp = malloc(n * sizeof(double));

    for (kick = 0; kick <= nkick; ++kick) {
        // Accumulate the data that are in the window before this kick.
        count = 0;
        start = ind;
        exitflag = 0;
        while (ind < n && !exitflag) {
            if (kick == nkick || t[ind] <= kicks[kick])
                ttmp[count++] = t[ind++];
            else exitflag = 1;
        }

        // Update the orbital parameters to the values that they should have
        // after this kick.
        for (i = 0; i < np; ++i) {
            e[i]      += del_e[kick*np+i];
            a[i]      += del_a[kick*np+i];
            t0[i]     += del_t0[kick*np+i];
            pomega[i] += del_pomega[kick*np+i];
            ix[i]     += del_ix[kick*np+i];
            iy[i]     += del_iy[kick*np+i];
        }

        // Solve the system using the Kepler solver.
        info = kepler_solve (count, ttmp, mstar, np, m, a, t0, e, pomega,
                             ix, iy, coords + (np+1)*6*start);

        // Set the orbital parameters back to their previous values.
        for (i = 0; i < np; ++i) {
            e[i]      -= del_e[kick*np+i];
            a[i]      -= del_a[kick*np+i];
            t0[i]     -= del_t0[kick*np+i];
            pomega[i] -= del_pomega[kick*np+i];
            ix[i]     -= del_ix[kick*np+i];
            iy[i]     -= del_iy[kick*np+i];
        }

        // Exit if the solve didn't succeed.
        if (info != 0) break;
    }

    free(ttmp);
    return info;
}
