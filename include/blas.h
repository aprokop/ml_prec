#ifndef _BLAS_H_
#define _BLAS_H_

#include "tools.h"

extern "C" {
	void	FORTRAN(dscal)  (const int* N, const double* ALPHA, double X[], const int* INCX);
	void	FORTRAN(dcopy)  (const int* N, const double X[], const int* INCX, double Y[], const int* INCY);
	void	FORTRAN(daxpy)  (const int* N, const double* ALPHA, const double X[], const int* INCX,
		double Y[], const int* INCY);
	double	FORTRAN(ddot)   (const int* N, const double X[], const int* INCX, const double Y[], const int* INCY);
	double	FORTRAN(dnrm2)  (const int* N, const double X[], const int* INCX);
	double	FORTRAN(dasum)  (const int* N, const double X[], const int* INCX);
	double  FORTRAN(idamax) (const int* N, const double X[], const int* INCX);
}

#endif // _BLAS_H_
