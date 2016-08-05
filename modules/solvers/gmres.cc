#include "solvers.h"
#include "include/logger.h"
#include "include/time.h"

#include <cmath>

DEFINE_LOGGER("GMRES");

static
void update_x(uint i, const PrecBase& B, const std::vector<Vector>& v, double** hh, uvector<double>& rs,
              Vector& x, Vector& tmp);

void GMRESSolver    (const CSRMatrix& A, const Vector& b, const PrecBase& B, Vector& x, SolverStats& stats,
                     double eps, NormType norm_type, bool silent) THROW {
    double  gtime = pclock();
    ASSERT_SIZE(b.size(), A.size());
    ASSERT_SIZE(x.size(), A.size());

    uint n = A.size();

    const uint m = 30;		    /* restart */

    uvector<double> s(m), c(m);	    /* sin and cos for rotations */
    uvector<double> rs(m);		    /* rhs for the system (implicit residual) */

    /* Hessenberg system related stuff */
    uvector<double> H((m+1)*m);
    double *hh[m];
    for (uint k = 0; k < m; k++)
        hh[k] = &H[k*m];

    /* Arnoldi's basis vectors */
    std::vector<Vector> v(m);
    for (uint k = 0; k < m; k++)
        v[k].resize(n);

    Vector z(n);

    double norm, init_norm;
    double dtmp;

    double  mult = 0,  inv = 0,  cstr = 0, delta;
    int	   nmult = 0, ninv = 0;

    generate_x0(x);
    residual(A, b, x, v[0]);
    norm = init_norm = dnrm2(v[0]);

#ifdef ABSOLUTE_NORM
    init_norm = 1;
#else
    check_and_replace_eps(init_norm, eps);
#endif

    int niter = 0;
    while (true) {
        if (norm/init_norm <= eps)
            return;

        dscal(1./norm, v[0]);

        /* Initialize 1st rhs term of H system */
        rs[0] = norm;

        uint i = 0;
        for (; i < m-1;) {
            /* v[i+1] = A M^-1 v[i] */
            delta = pclock();
            B.solve(v[i], z);
            inv += pclock() - delta;
            ninv++;

            delta = pclock();
            multiply(A, z, v[i+1]);
            mult += pclock() - delta;
            nmult++;

            /* Orhogonalize */
            if (true) {
                /* Classical Gram-Schmidt */
                for (uint k = 0; k <= i; k++)
                    hh[k][i] = ddot(v[i+1], v[k]);
                for (uint k = 0; k <= i; k++)
                    daxpy(-hh[k][i], v[k], v[i+1]);
            } else {
                /* Modified Gram-Schmidt */
                for (uint k = 0; k <= i; k++) {
                    hh[k][i] = ddot(v[i+1], v[k]);
                    daxpy(-hh[k][i], v[k], v[i+1]);
                }
            }

            /* Normalize */
            hh[i+1][i] = dnrm2(v[i+1]);
            dscal(1./hh[i+1][i], v[i+1]);

            /* Update factorization of hh by plane rotation */
            for (uint k = 1; k <= i; k++) {
                int k1 = k-1;
                dtmp = hh[k1][i];
                hh[k1][i] =  c[k1]*dtmp + s[k1]*hh[k][i];
                hh[k][i]  = -s[k1]*dtmp + c[k1]*hh[k][i];
            }

            /* Determine next plane rotation */
            dtmp = sqrt(hh[i][i]*hh[i][i] + hh[i+1][i] * hh[i+1][i]);
            dtmp = 1.0/dtmp;

            c[i]    =  hh[i][i]   * dtmp;
            s[i]    =  hh[i+1][i] * dtmp;

            /* Perform final rotation */
            hh[i][i] = c[i]*hh[i][i] + s[i]*hh[i+1][i];

            /* Update rhs */
            rs[i+1]  = -s[i]*rs[i];
            rs[i]   *=  c[i];

            /* Determine residual norm and test convergence */
            norm = fabs(rs[i+1]);

            niter++;
            i++;

            if (!silent)
                LOG_DEBUG("#" << niter << ": relative -> " << std::scientific << norm/init_norm << "   absolute -> " << norm);

            if (norm/init_norm <= eps)
                break;
        }

        /* Update solution vector */
        update_x(i, B, v, hh, rs, x, z);

        if (norm/init_norm <= eps)
            return;

        residual(A, b, x, v[0]);
    }

}

/* n - dimension of the current Krylov space */
void update_x(uint n, const PrecBase& B, const std::vector<Vector>& v, double** hh, uvector<double>& rs,
              Vector& x, Vector& z) {
    if (n == 0)
        return;
    n--;

    /* Solve upper triangular system */

    rs[n] /= hh[n][n];
    for (uint i = 1; i <= n; i++) {
        int k = n-i;	/* row index */

        for (uint j = k+1; j < n; j++)
            rs[k] -= hh[k][j] * rs[j];
        rs[k] /= hh[k][k];
    }

    /* Form linear combination to get solution */
    Vector tmp(x.size(), 0.0);
    for (uint j = 0; j <= n; j++)
        daxpy(rs[j], v[j], tmp);
    B.solve(tmp, z);

    daxpy(1.0, z, x);
}

