#ifndef _NOBODY_BLAS_H_
#define _NOBODY_BLAS_H_

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define NBODY_SOL2

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define NBODY_SGI

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define NBODY_LINUX

#elif defined (__APPLE__)
#define NBODY_MAC

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define NBODY_AIX

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define NBODY_ALPHA

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define NBODY_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define NBODY_CYGWIN
#else
#define NBODY_WINDOWS
#define BLAS_NO_UNDERSCORE
#endif

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || \
      defined (ARCH_HPUX) || defined (__hp700) || defined (MHP700) || \
      defined (ARCH_HP700)
#define NBODY_HP
#define BLAS_NO_UNDERSCORE

#endif

#if defined (BLAS_NO_UNDERSCORE)

#define BLAS_DSWAP dswap
#define BLAS_DCOPY dcopy
#define BLAS_DAXPY daxpy
#define BLAS_DSCAL dscal
#define BLAS_DNRM2 dnrm2
#define BLAS_IDAMAX idamax

#else

#define BLAS_DSWAP dswap_
#define BLAS_DCOPY dcopy_
#define BLAS_DAXPY daxpy_
#define BLAS_DSCAL dscal_
#define BLAS_DNRM2 dnrm2_
#define BLAS_IDAMAX idamax_

#endif

void BLAS_DSWAP (int *n, double *dx, int *incx, double *dy, int *incy);
#define BLAS_dswap(n,dx,incx,dy,incy) \
{ \
    int N = n, INCX = incx, INCY = incy; \
    BLAS_DSWAP (&N, dx, &INCX, dy, &INCY); \
}

void BLAS_DCOPY (int *n, double *dx, int *incx, double *dy, int *incy);
#define BLAS_dcopy(n,dx,incx,dy,incy) \
{ \
    int N = n, INCX = incx, INCY = incy; \
    BLAS_DCOPY (&N, dx, &INCX, dy, &INCY); \
}

void BLAS_DAXPY (int *n, double *da, double *dx, int *incx, double *dy,
                 int *incy);
#define BLAS_daxpy(n,da,dx,incx,dy,incy) \
{ \
    double DA = da; \
    int N = n, INCX = incx, INCY = incy; \
    BLAS_DAXPY (&N, &DA, dx, &INCX, dy, &INCY); \
}

void BLAS_DSCAL (int *n, double *da, double *dx, int *incx);
#define BLAS_dscal(n,da,dx,incx) \
{ \
    double DA = da; \
    int N = n, INCX = incx; \
    BLAS_DSCAL (&N, &DA, dx, &INCX); \
}

double BLAS_DNRM2 (int *n, double *dx, int *incx);
int BLAS_IDAMAX (int *n, double *dx, int *incx);

#endif
