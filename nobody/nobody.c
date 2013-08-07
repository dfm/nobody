#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int nparticles;
double *masses;

void setup (int n, double *m)
{
    int i;
    nparticles = n;
    masses = malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) masses[i] = m[i];
}

void cleanup ()
{
    free(masses);
}

void gradient (double t, double *x, double *xout)
{
    int i, j, k;
    double delta[3], d2;

    for (i = 0; i < nparticles; ++i) {
        xout[6*i]   = x[6*i+3];
        xout[6*i+1] = x[6*i+4];
        xout[6*i+2] = x[6*i+5];
        xout[6*i+3] = 0.0;
        xout[6*i+4] = 0.0;
        xout[6*i+5] = 0.0;
    }

    for (i = 0; i < nparticles; ++i) {
        for (j = i + 1; j < nparticles; ++j) {
            d2 = 0.0;
            for (k = 0; k < 3; ++k) {
                delta[k] = x[6*i+k] - x[6*j+k];
                d2 += delta[k] * delta[k];
            }
            for (k = 0; k < 3; ++k) {
                delta[k] /= d2 * sqrt(d2);
                xout[6*i+3+k] -= masses[j] * delta[k];
                xout[6*j+3+k] += masses[i] * delta[k];
            }
        }
    }
}

double midpoint_step (int n, double t, double *x, double H, int nstep,
                      double *xout,
                      void (*dxdt) (double, double*, double*))
{
    int i, k;
    double h = H / nstep, swap,
           *x0 = malloc(n * sizeof(double)),
           *x1 = malloc(n * sizeof(double)),
           *dx = malloc(n * sizeof(double));

    dxdt(t, x, dx);
    for (i = 0; i < n; ++i) {
        x0[i] = x[i];
        x1[i] = x[i] + h * dx[i];
    }

    for (k = 2; k <= nstep; ++k) {
        t += h;
        dxdt(t, x1, dx);
        for (i = 0; i < n; ++i) {
            swap = x1[i];
            x1[i] = x0[i] + 2 * h * dx[i];
            x0[i] = swap;
        }
    }

    t += h;
    dxdt(t, x1, dx);
    for (i = 0; i < n; ++i)
        xout[i] = 0.5 * (x0[i] + x1[i] + h * dx[i]);

    free(x0);
    free(x1);
    free(dx);

    return t;
}

void polynomial_extrapolate_zero (int ind, int n, double *h, double *x,
                                  double **tableau, double *xout, double *xerr)
{
    int i, k;
    double *c, delta, f1, f2;

    for (i = 0; i < n; ++i) xout[i] = xerr[i] = x[i];
    if (ind == 0) {
        for (i = 0; i < n; ++i) tableau[ind][i] = x[i];
        return;
    }

    c = malloc(n * sizeof(double));
    for (i = 0; i < n; ++i) c[i] = x[i];
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

    for (i = 0; i < n; ++i) tableau[ind][i] = xerr[i];
    free(c);
}

double bs_step (int n, double t, double *x, double H, double eps,
                double *xfinal,
                void (*dxdt) (double, double*, double*))
{
    int i, k, kmax = 8, nstep;
    double err, maxerr, t0 = t,
           *h = malloc(kmax * sizeof(double)),
           *xest = malloc(n * sizeof(double)),
           *xout = malloc(n * sizeof(double)),
           *xerr = malloc(n * sizeof(double)),
           **tableau = malloc(kmax * sizeof(double*));
    for (k = 0; k < kmax; ++k) {
        nstep = 2 * (k + 1);

        t = midpoint_step (n, t0, x, H, nstep, xest, dxdt);
        h[k] = H * H / nstep / nstep;
        tableau[k] = malloc(n * sizeof(double));
        polynomial_extrapolate_zero (k, n, h, xest, tableau, xout, xerr);

        maxerr = xerr[0] * xerr[0];
        for (i = 1; i < n; ++i) {
            err = xerr[i] * xerr[i];
            if (err > maxerr) maxerr = err;
        }
        maxerr = sqrt(maxerr);
        if (maxerr <= eps) break;
    }

    free(h);
    free(xest);
    free(xerr);
    for (i = 0; i < k; ++i) free(tableau[i]);
    free(tableau);

    if (maxerr > eps) {
        if (0.5 * H > 0.0) {
            free(xout);
            return bs_step(n, t0, x, 0.5 * H, eps, xfinal, dxdt);
        }
        fprintf(stderr, "Underflow encountered in BS step size.");
    }

    for (i = 0; i < n; ++i) xfinal[i] = xout[i];
    free(xout);
    return t;
}

int main ()
{
    int i, j, k;
    double m[2]  = {1.0, 2.0},
           x[12] = {1.0, 0.0, 0.0, 0.0,  0.5, 0.0,
                    0.0, 0.0, 0.0, 0.0, -0.1, 0.0},
           g[12], h = 0.5, t = 0.0;

    setup(2, m);

    for (k = 0; k < 1000; ++k) {
        t = bs_step (12, t, x, h, 1e-9, x, &gradient);
        printf("%e ", t);
        for (i = 0; i < 12; ++i) printf("%e ", x[i]);
        printf("\n");
    }

    return 0;
}
