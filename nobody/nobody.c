#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define KMAX 8

//
// Global definitions for the N-body system.
//
int nparticles;
double *masses;

//
// Set up the system. This saves the number of particles and the particle
// masses to the globals defined above for the computation of the
// accelerations.
//
// :param int n:       The number of particles in the system.
// :param double m[n]: The masses of the particles.
//
void nobody_setup (int n, double *m)
{
    int i;
    nparticles = n;
    masses = malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) masses[i] = m[i];
}

//
// Clean up the memory usage of the system.
//
void nobody_cleanup ()
{
    free(masses);
}

//
// Compute the phase space gradient of a given N-body system using direct
// summation. The masses are stored in the global list ``masses`` but the
// phase space coordinates should be provided in the form
// ``(x_1, y_1, z_1, v_{x,1}, v_{y,1}, v_{z_1}, ...)`` and the computed
// gradient will be returned in the form
// ``(v_{x,1}, v_{y,1}, v_{z_1}, a_{x,1}, a_{y,1}, a_{z_1}, ...)``.
//
// :param double t:                  The current time step (unused).
// :param double x[6*nparticles]:    The current phase space coordinates.
// :param double xout[6*nparticles]: The computed accelerations. (output)
//
void nobody_gradient (double t, double *x, double *xout)
{
    int i, j, k;
    double delta[3], d2;

    // Copy the velocities and initialize the accelerations.
    for (i = 0; i < nparticles; ++i) {
        xout[6*i]   = x[6*i+3];
        xout[6*i+1] = x[6*i+4];
        xout[6*i+2] = x[6*i+5];
        xout[6*i+3] = 0.0;
        xout[6*i+4] = 0.0;
        xout[6*i+5] = 0.0;
    }

    // Run a direct loop over pairs of particles.
    for (i = 0; i < nparticles; ++i) {
        for (j = i + 1; j < nparticles; ++j) {
            // Compute the separation vector and distance between the
            // particles.
            d2 = 0.0;
            for (k = 0; k < 3; ++k) {
                delta[k] = x[6*i+k] - x[6*j+k];
                d2 += delta[k] * delta[k];
            }

            // Update the accelerations based on the reflex gravitational
            // force.
            for (k = 0; k < 3; ++k) {
                delta[k] /= d2 * sqrt(d2);
                xout[6*i+3+k] -= masses[j] * delta[k];
                xout[6*j+3+k] += masses[i] * delta[k];
            }
        }
    }
}

//
// Integrate a first order ODE using the midpoint method.
//
// :param int n:          The dimension of the input space.
// :param double t:       The current time step.
// :param double x[n]:    The initial coordinates.
// :param double H:       The amount of time to advance the integration for.
// :param int nstep:      The number of steps to use to cross ``H``.
// :param double xout[n]: The final coordinates. (output)
// :param void (*dxdt):   A function to compute the gradient of ``x`` at a
//                        particular value of ``x``.
//
// :returns double t:     The final time of the system.
//
double nobody_midpoint (int n, double t, double *x, double H, int nstep,
                        double *xout,
                        void (*dxdt) (double, double*, double*))
{
    int i, k;
    double h = H / nstep, swap,
           *x0 = malloc(n * sizeof(double)),
           *x1 = malloc(n * sizeof(double)),
           *dx = malloc(n * sizeof(double));

    // Compute the initial gradient and save the first two terms in the
    // stencil using forward Euler.
    dxdt(t, x, dx);
    for (i = 0; i < n; ++i) {
        x0[i] = x[i];
        x1[i] = x[i] + h * dx[i];
    }

    // Loop over the requested number of midpoint steps and propagate the
    // gradient using centered difference.
    for (k = 2; k <= nstep; ++k) {
        t += h;
        dxdt(t, x1, dx);
        for (i = 0; i < n; ++i) {
            swap = x1[i];
            x1[i] = x0[i] + 2 * h * dx[i];
            x0[i] = swap;
        }
    }

    // Compute the final position based on the last centered difference.
    t += h;
    dxdt(t, x1, dx);
    for (i = 0; i < n; ++i)
        xout[i] = 0.5 * (x0[i] + x1[i] + h * dx[i]);

    // Clean up.
    free(x0);
    free(x1);
    free(dx);

    return t;
}

//
// Use polynomial extrapolation to extrapolate a function to zero. This is
// based on the implementation of Neville's algorithm given in Numerical
// Recipes. This should be called iteratively for decreasing values of the
// independent variable (h) and the estimate will be compounded in the array
// called tableau. This function also computes a conservative error estimate
// on the extrapolated value.
//
// :param int ind:         The iteration number.
// :param int n:           The number of functions that you're trying to
//                         extrapolate.
// :param double h[>=ind]: An array containing the current value of the
//                         independent variable in the ``ind``-th entry and
//                         all previous (larger) values.
// :param double x[n]:     The estimated value of the functions at the current
//                         value of ``h``.
// :param double tableau[>=ind][n]:
//                         The tableau of coefficients that have been computed
//                         so far.
// :param double xout[n]:  The current estimate of the functions evaluated at
//                         ``h = 0``.
// :param double xerr[n]:  A conservative estimate of the error in ``xout``.
//                         This value might be negative so you should take the
//                         absolute value or the square before interpreting.
//
void nobody_poly_extrap_zero (int ind, int n, double *h, double *x,
                              double **tableau, double *xout, double *xerr)
{
    int i, k;
    double *c, delta, f1, f2;

    // Initialize the output at the current value.
    for (i = 0; i < n; ++i) xout[i] = xerr[i] = x[i];

    // On the first pass, save the initial estimate directly to the tableau.
    if (ind == 0) {
        for (i = 0; i < n; ++i) tableau[ind][i] = x[i];
        return;
    }

    c = malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) c[i] = x[i];

    // Traverse the tableau and reduce the output estimate.
    for (k = 1; k < ind; ++k) {
        delta = 1.0 / (h[ind-k] - h[ind]);
        f1 = h[ind] * delta;
        f2 = h[ind-k] * delta;
        for (i = 0; i < n; ++i) {
            delta = c[i] - tableau[k][i];
            tableau[k][i] = xerr[i];
            xerr[i] = f1 * delta;
            c[i] = f2 * delta;
            xout[i] += xerr[i];
        }
    }

    // Update the current column in the tableau.
    for (i = 0; i < n; ++i) tableau[ind][i] = xerr[i];
    free(c);
}

//
// Integrate a time step in a first order ODE using the Bulirsch–Stoer
// algorithm. Briefly, the algorithm works by iteratively running the midpoint
// method with increasing numbers of midpoint steps and then extrapolating the
// result to an infinite number of midpoint steps. If the required accuracy
// is not reached by using (at most) 16 midpoint steps, a smaller time step
// is tried iteratively until convergence is reached. This implementation
// isn't particularly intelligent: it simply tries halving the time step until
// it converges of reaches a time step that is zero to machine precision.
// This function returns the final time of the system after applying the final
// converged time step.
//
// :param int n:          The dimension of the input space.
// :param double t:       The current time step.
// :param double x[n]:    The initial coordinates.
// :param double H:       The time step to try. There is no guarantee that
//                        the system will *actually* be integrated for this
//                        amount of time.
// :param double eps:     The required extrapolation accuracy.
// :param double xout[n]: The final coordinates. (output)
// :param void (*dxdt):   A function to compute the gradient of ``x`` at a
//                        particular value of ``x``.
//
// :returns double t:     The final time of the system. In general, this
//                        *won't* be equal to the initial time plus ``H``
//                        because a smaller time step might be required to
//                        achieve the target accuracy ``eps``.
//
double nobody_bs (int n, double t, double *x, double H, double eps,
                  double *xfinal, void (*dxdt) (double, double*, double*))
{
    int i, k, kmax = 8, nstep;
    double err, maxerr, t0 = t, h[KMAX], *tableau[KMAX],
           *xest = malloc(n * sizeof(double)),
           *xout = malloc(n * sizeof(double)),
           *xerr = malloc(n * sizeof(double));

    // Loop over possible numbers of midpoint steps.
    for (k = 0; k < KMAX; ++k) {
        // Run a midpoint integration for ``nstep = 2*(k+1)`` steps.
        nstep = 2 * (k + 1);
        t = nobody_midpoint (n, t0, x, H, nstep, xest, dxdt);

        // Attempt a polynomial extrapolation to zero step size.
        h[k] = H * H / nstep / nstep;
        tableau[k] = malloc(n * sizeof(double));
        nobody_poly_extrap_zero (k, n, h, xest, tableau, xout, xerr);

        // Find the maximum error estimate in the system.
        maxerr = xerr[0] * xerr[0];
        for (i = 1; i < n; ++i) {
            err = xerr[i] * xerr[i];
            if (err > maxerr) maxerr = err;
        }

        // Check for convergence.
        maxerr = sqrt(maxerr);
        if (maxerr <= eps) break;
    }

    // Clean up.
    free(xest);
    free(xerr);
    for (i = 0; i < k; ++i) free(tableau[i]);

    // If the system did not converge to the required accuracy, try a smaller
    // time step.
    if (maxerr > eps) {
        t = nobody_bs (n, t0, x, 0.5 * H, eps, xfinal, dxdt);
        if (t > t0) {
            free(xout);
            return t;
        }

        // If we get here, the step size was zero to machine precision. We'll
        // warn and return the most recent attempt at the integration.
        fprintf(stderr, "Underflow encountered in BS step size.\n");
    }

    // Copy the results of the integration to the output array.
    for (i = 0; i < n; ++i) xfinal[i] = xout[i];
    free(xout);

    // Return the final time of the system.
    return t;
}

// int main ()
// {
//     int i, j, k;
//     double m[2]  = {1.0, 2.0},
//            x[12] = {1.0, 0.0, 0.0, 0.0,  0.5, 0.0,
//                     0.0, 0.0, 0.0, 0.0, -0.1, 0.0},
//            g[12], h = 0.5, t = 0.0;
//
//     nobody_setup(2, m);
//
//     for (k = 0; k < 1000; ++k) {
//         t = nobody_bs (12, t, x, h, 1e-9, x, &nobody_gradient);
//         printf("%e ", t);
//         for (i = 0; i < 12; ++i) printf("%e ", x[i]);
//         printf("\n");
//     }
//
//     return 0;
// }
