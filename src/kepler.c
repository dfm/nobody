#include "nobody.h"

//
// The maximum number of iterations to attempt when trying to solve Kepler's
// equation.
//
#define MAXITER 500

//
// A MAGIC number setting the convergence tolerance for Kepler's equation.
//
#define CONVTOL 1.48e-10

//
// A tiny number that is approximately machine precision.
//
#define EPSILON 2e-16

//
// Get the mean anomaly associated with a particular eccentric anomaly and
// eccentricity.
//
// :param double psi:   The eccentric anomaly.
// :param double e:     The eccentricity.
//
// :returns double psi: The mean anomaly.
//
double kepler_ecc_to_mean_anomaly (double psi, double e)
{
    return psi - e * sin(psi);
}

//
// Iteratively solve Kepler's equation to find the eccentric anomaly of an
// orbit given a mean anomaly and an eccentricity.
//
// :param double M:     The mean anomaly.
// :param double e:     The eccentricity.
// :param int *info:    The output success flag. 0=success, 1=did not
//                      converge after 500 iterations, 2=unphysical
//                      eccentricity.
//
// :returns double psi: The eccentric anomaly.
//
double kepler_mean_to_ecc_anomaly (double M, double e, int *info)
{
    int i;
    double wt0, psi0, psi = 0.0;

    // Check for un-physical parameters.
    if (e < 0 || e >= 1) {
        *info = 2;
        return 0.0;
    }

    *info = 0;
    wt0 = fmod(M, 2 * M_PI);
    psi0 = wt0;

    for (i = 0; i < MAXITER; ++i) {
        // Compute the function value and the first two derivatives.
        double spsi0 = sin(psi0),
               f = psi0 - e * spsi0 - wt0,
               fp = 1.0 - e * cos(psi0),
               fpp = e * spsi0;

        // Take a second order step.
        psi = psi0 - 2.0 * f * fp / (2.0 * fp * fp - f * fpp);

        // Accept?
        if (fabs(psi - psi0) <= CONVTOL) return psi;

        // Save as the previous step.
        psi0 = psi;
    }

    // If we get here, we didn't ever converge.
    *info = 1;
    return psi;
}

//
// Solve Kepler's equation for the orbit of one planet.
//
// :param int n:              The number of time steps to solve.
// :param double t[n]:        The times where the orbit should be computed.
// :param double mstar:       The mass of the star in Solar masses.
// :param double mplanet:     The mass of the planet in Solar masses.
// :param double a:           The semi-major axis of the orbit in Solar radii.
// :param double t0:          The epoch of the orbit (defined as the center
//                            time of a target transit) in days.
// :param double e:           The eccentricity of the orbit.
// :param double pomega:      The orientation of the orbital ellipse. Defined
//                            as the angle between the major-axis of the
//                            ellipse and the observer in radians.
// :param double ix:          The inclination of the orbit relative to the
//                            observer in radians. At ``ix=0``, the system is
//                            edge on and at ``ix=0.5*M_PI``, the system is
//                            face on.
// :param double iy:          The orientation of the orbit in the plane of the
//                            sky measured in radians.
// :param double coords[6*n]: The phase space coordinates of the planet at
//                            each time ``t``. The coordinate system is such
//                            that the first element points towards the
//                            observer and the other two are in the plane of
//                            the sky (as defined by ``iy``). The positions
//                            are measured in Solar radii and the velocities
//                            in Solar radii per day.
//
// :returns int info:         A flag describing the success of the computation.
//                            0=success, 1=solve did not converge after 500
//                            iterations, 2=unphysical eccentricity.
//
int kepler_solve_one (int n, double *t, double mstar, double mplanet,
                      double a, double t0, double e, double pomega, double ix,
                      double iy, double *xout, double *vout, int incr)
{
    int i, offset, info = 0, sgn;
    double period = 2 * M_PI * sqrt(a * a * a / G_GRAV / (mstar + mplanet)),
           psi0 = 2 * atan2(tan(0.5 * pomega), sqrt((1 + e) / (1 - e))),
           t1 = t0 -  0.5 * period * (psi0 - e * sin(psi0)) / M_PI,
           manom, psi, cpsi, spsi, d, cth, sth, r, x, y, xp, yp, xsx,
           dmanomdt, denom, drdt, dthdt, dxdt, dydt, dxpdt, dypdt,
           dxsxdt, cpomega = cos(pomega), spomega = sin(pomega),
           cix = cos(ix), six = sin(ix), ciy = cos(iy), siy = sin(iy);

    if (e < 0.0 || fabs(e - 1.0) <= EPSILON) return 2;

    dmanomdt = 2 * M_PI / period;
    for (i = 0; i < n; ++i) {
        manom = dmanomdt * (t[i] - t1);
        psi = kepler_mean_to_ecc_anomaly (manom, e, &info);

        // Did solve fail?
        if (info != 0) return info;

        // Compute the derivative of the eccentric anomaly.
        cpsi = cos(psi);
        spsi = sin(psi);

        // Compute the true anomaly.
        d = 1.0 - e * cpsi;
        if (d == 0.0) return 3;
        cth = (cpsi - e) / d;
        sth = sqrt(1 - cth * cth);

        // Compute the radius and derivative.
        denom = 1.0 / (1 + e * cth);
        r = a * (1 - e * e) * denom;
        denom = sqrt(1 - e * e);
        drdt = dmanomdt * a * e * sth / denom;
        dthdt = dmanomdt * a * a * denom / r / r;

        // Compute the coordinates in the plane of the orbit.
        sgn = (spsi >= 0) - (spsi < 0);
        x = r * cth;
        y = r * sth * sgn;
        dxdt = drdt * cth - r * sth * dthdt;
        dydt = (drdt * sth + r * cth * dthdt) * sgn;

        // Rotate by pomega.
        xp =  x * cpomega + y * spomega;
        yp = -x * spomega + y * cpomega;

        // Rotate by the inclination angle.
        xsx = xp * six;

        // Compute the positions.
        offset = incr*i;
        xout[offset]   = xp * cix;
        xout[offset+1] = yp * ciy - xsx * siy;
        xout[offset+2] = yp * siy + xsx * ciy;

        // Rotate the velocities.
        dxpdt =  dxdt * cpomega + dydt * spomega;
        dypdt = -dxdt * spomega + dydt * cpomega;
        dxsxdt = dxpdt * six;

        // Compute the velocities.
        vout[offset]   = sgn * (dxpdt * cix);
        vout[offset+1] = sgn * (dypdt * ciy - dxsxdt * siy);
        vout[offset+2] = sgn * (dypdt * siy + dxsxdt * ciy);
    }

    return info;
}

//
// Solve Kepler's equation for the orbits of a system of planets.
//
// :param int n:             The number of time steps to solve.
// :param double t[n]:       The times where the orbit should be computed.
// :param double mstar:      The mass of the star in Solar masses.
// :param int np:            The number of planets in the system.
// :param double m[np]:      The masses of the planets in Solar masses.
// :param double a[np]:      The semi-major axis of the orbit in Solar radii.
// :param double t0[np]:     The epochs of the orbits (defined as the center
//                           time of a target transit) in days.
// :param double e[np]:      The eccentricities of the orbits.
// :param double pomega[np]: The orientations of the orbital ellipses. Defined
//                           as the angle between the major-axis of the
//                           ellipse and the observer in radians.
// :param double ix[np]:     The inclinations of the orbits relative to the
//                           observer in radians. At ``ix=0``, the system is
//                           edge on and at ``ix=0.5*M_PI``, the system is
//                           face on.
// :param double iy[np]:     The orientations of the orbits in the plane of
//                           the sky measured in radians.
// :param double coords[6*(np+1)*n]:
//                        The phase space coordinates of the star and the
//                        planets at each time ``t``. The coordinate system
//                        is such that the first element points towards the
//                        observer and the other two are in the plane of
//                        the sky (as defined by ``iy``). The positions
//                        are measured in Solar radii and the velocities
//                        in Solar radii per day.
//
// :returns int info:     A flag describing the success of the computation.
//                        0=success, 1=solve did not converge after 500
//                        iterations, 2=unphysical eccentricity.
//
int kepler_solve (int n, double *t, double mstar, int np, double *m,
                  double *a, double *t0, double *e, double *pomega,
                  double *ix, double *iy, double *x, double *v)
{
    int i, j, k, offset, info = 0;

    // Initialize the star's coordinates;
    offset = 3 * (np + 1);
    for (i = 0; i < n; ++i)
        for (j = 0; j < 3; ++j) {
            x[offset*i+j] = 0.0;
            v[offset*i+j] = 0.0;
        }

    // Loop over the planets and solve.
    for (k = 0; k < np; ++k) {
        // Solve the system for one planet.
        info = kepler_solve_one (n, t, mstar, m[k], a[k], t0[k], e[k],
                                 pomega[k], ix[k], iy[k],
                                 x+3*(1+k), v+3*(1+k), offset);
        if (info != 0) {
            return info;
        }

        // Copy over the results.
        for (i = 0; i < n; ++i) {
            // Update the stellar velocity.
            for (j = 0; j < 3; ++j)
                v[offset*i+j] -= m[k] * v[offset*i+j+3*(k+1)] / mstar;
        }
    }

    return info;
}
