#include "include/logger.h"
#include "include/fortran.h"
#include "amg_prec.h"

DEFINE_LOGGER("AMGPrec");

extern "C" {
    void FORTRAN(amg1r6)(double *a, int *ia, int *ja,  double *u, double *f, int *ig,
                         int *nda, int *ndia, int *ndja, int *ndu, int *ndf, int *ndig,
                         int *nnu, int *matrix, int *iswtch, int *iout, int *iprint,
                         int *levelx, int *ifirst, int *ncyc, double *eps, int *madapt, int *nrd, int *nsolco,
                         int *nru, double *ecg1, double *ecg2, double *ewt2, int *nwt, int *ntr, int *ierr);
}

void AMGPrec::solve(Vector& f, Vector& x) const THROW {
    ASSERT((int)f.size() == n && (int)x.size() == n, "Wrong dimension: n = " << n << ", f = " << f.size() << ", x = " << x.size());

    x.resize(ndu);
    f.resize(ndf);

    int ierr;
    FORTRAN(amg1r6)(&a[0], &ia[0], &ja[0], &x[0], &f[0], &ig[0],
                    &nda, &ndia, &ndja, &ndu, &ndf, &ndig, &n,
                    &amg_config.matrix,
                    &amg_config.iswtch,
                    &amg_config.iout,
                    &amg_config.iprint,
                    &amg_config.levelx,
                    &amg_config.ifirst,
                    &amg_config.ncyc,
                    &amg_config.eps,
                    &amg_config.madapt,
                    &amg_config.nrd,
                    &amg_config.nsolco,
                    &amg_config.nru,
                    &amg_config.ecg1,
                    &amg_config.ecg2,
                    &amg_config.ewt2,
                    &amg_config.nwt,
                    &amg_config.ntr,
                    &ierr);

    // no setup after the first call
    amg_config.iswtch = 3;

    x.resize(n);
    f.resize(n);
}

AMGPrec::AMGPrec(const SkylineMatrix& A) {
    n = A.rows();
    ASSERT(A.cols() == A.rows(), "Matrix must be square");
    int ind = 0;

    ia.resize(A.ia.size());
    for (uint i = 0; i < A.ia.size(); i++) {
        ia[i] = A.ia[i] + 1;
    }
    ja.resize(A.ja.size());
    for (uint j = 0; j < A.ja.size(); j++)
        ja[j] = A.ja[j] + 1;
    a = A.a;

    nnz = A.nnz();
    LOG_DEBUG("nnz = " << nnz);

#if 0
    ndia = int(2.2*n);
    ndja = 3*n + 5*nnz;
    nda  = 3*n + 5*nnz;
    ndig = int(5.4*n);
    ndu  = int(2.2*nnz);
    ndf  = int(2.2*nnz);
#else
    // these values are sufficient for SPE test
    ndia =  5000000;
    ndja = 50000000;
    nda  = 50000000;
    ndig =  7000000;
    ndu  =  5000000;
    ndf  =  5000000;
#endif

    ia.resize(ndia);
    ja.resize(ndja);
    a.resize (nda);
    ig.resize(ndig);
}

AMGPrec::AMGConfig::AMGConfig() {
    /******************************************************************
     *                            GROUP 1
     ******************************************************************/
    /**
     * matrix - integer value containing info about matrix
     * 1st digit:
     *	    1 - symmetric,
     *	    2 - not symmetric
     * 2nd digit:
     *	    1 - rowsum zero,
     *	    2 - rowsum not zero
     */
    // matrix = 12;
    matrix = 22;

    /******************************************************************
     *                            GROUP 2
     ******************************************************************/
    /**
     * iswtch - parameter controlling which modules of amg1r6 are to be used
     *	    1 - -----, -----, -----, wrkcnt
     *	    2 - -----, -----, solve, wrkcnt
     *	    3 - -----, first, solve, wrkcnt
     *	    4 - setup, first, solve, wrkcnt
     *	    where
     *		setup  - defines the operators needed in the solution phase
     *		first  - initializes the solution vector (see ifirst)
     *		solve  - computes the solution by AMG cycling (see ncyc)
     *		wrkcnt - provides the user with info aout residuals, storage
     *			 requirements and CP-times (see iout)
     *	    if amg1r6 is called the first time, iswtch has to be 4.
     */
    iswtch = 4;

    /**
     * iout - parameter controlling the amount of output during solution phase
     * 1st digit:   not used, has to be non-zero
     * 2nd digit:
     *	    0 - no output,
     *	    1 - residual before and after solution process,
     *	    2 - add. statistics on CP-times and storage requirements
     *	    3 - add. residual after each AMG cycle
     */
    iout   = 10;

    /**
     * iprint - parameter specifying the FORTRAN unit numbers for output
     * 1st digit:   not used; has to be non-zero
     * 2nd & 3d :   unit number for results
     * 4th & 5th:   unit number for messages
     */
    iprint = 10606;	    // console output

    /******************************************************************
     *                            GROUP 3
     ******************************************************************/
    /**
     * maximum number of layers to be created (>=1)
     */
    levelx = 25;

    /**
     * ifirst - parameter for first approximation
     * 1st digit:	    not used; has to be non-zero
     * 2nd digit (itypu):
     *	    0 - no setting of first approximation
     *	    1 - first approximation constant to zero
     *	    2 - first approximation constant to one
     *	    3 - first approximation is random function with concrete random
     *		sequence being determined by the following digits
     * rest (rndu):
     *	    determines the concrete random sequence used in the case itypu=3
     *	    (ifirst=13 is equivalent to ifirst=1372815)
     */
    ifirst = 11;

    /**
     * ncyc - parameter describing the type of cycle to be used and the number
     *	      of cycles to be performed
     * 1st digit (igam):
     *	    1 - V-cycle
     *	    2 - V*-cycle
     *	    3 - F-cycle
     *	    4 - W-cycle
     *	    if ncyc is negative, then the approximation of the problem on the
     *	    second finest grid is computed by igam V-cycles on that particular
     *	    grid
     * 2nd digit (icgr):
     *	    0 - no Conjugate Gradient (CG)
     *	    1 - CG (only first step of CG)
     *	    2 - CG (full CG)
     * 3rd digit (iconv):
     *	    convergence criterion for the user-defined problem (finest grid):
     *	    1 - perform a fixed number of cycles as given by ncycle (see below)
     *	    2 - stop, if ||res|| < eps
     *	    3 - stop, if ||res|| < eps * |F|
     *	    4 - stop, if ||res|| < eps * |U| * |diag|,
     *	    with
     *		||REC|| - L2-norm of residual
     *		eps	- (see input parameter eps)
     *		|F|	- sup norm of rhs
     *		|U|	- sup norm of solution
     *		|diag|	- maximal diagonal entry in matrix
     *	    in any case the solution process stops after at most ncycle cycles
     * rest (ncycle):
     *	    maximal number of cycles to be performed (>0) or ncycle=0: no cycling
     */
    // ncyc = 122100;
    // ncyc = 102100;
    ncyc = 1011;

    /**
     * eps - convergence criterion for solution process (see parameter ncyc)
     * [Note: no more than NCYCLE cycles are performed, regardless of eps]
     */
    eps    = 1e-6;

    /**
     * madapt - integer value specifying the choice of coarsest grid in cycling
     * 1st digit (msel):
     *	    1 - in cycling, all grids constructed in the setup phase are used
     *		without check
     *	    2 - the number of grids is automatically reduced if the convergence
     *		factor on the coarser grids is found to be large than a given
     *		value fac (see below)
     * rest (fac):
     *	    the rest of madapt defines the fractional part of a real number fac
     *	    between 0.1 and 0.99; if madapt consists of only one digit, fac is
     *	    set to 0.7 by default
     */
    // madapt = 27;
    madapt = 17;

    /**
     * nrd - parameter describing relaxation (downwards)
     * 1st digit:	    not used; has to be non-zero
     * 2nd digit (nrdx):
     *	    actual number of smoothing steps to be performed
     *	    the type of which is given by the following digits
     * rest (nrdtyp):
     *	    1 - relaxation over the F-points only
     *	    2 - full GS sweep
     *	    3 - relaxation over the C-points only
     *	    4 - full more color sweep, highest color first
     */
    // nrd = 1131;
    nrd = 112;

    /**
     * nsolco - parameter controlling the solution on coarsest grid
     * 1st digit (nsc):
     *	    1 - Gauss-Seidel method
     *	    2 - direct solver (Yale smp)
     * rest (nrcx, only if nsc == 1):
     *	    number of GS sweeps on coarsest grid (>=0)
     *	    [if = 0 then as many as are needed to reduce residual
     *	     by two orders of magnitude]
     */
    nsolco = 2;

    /**
     * nru - parameter for relaxation (upwards), analogous to nrd
     */
    nru = nrd;

    /******************************************************************
     *                            GROUP 4
     ******************************************************************/
    /**
     * ecg1, ecg2, ewt2 - parameters affecting the creation of coarser grids and/or the
     * definition of the interpolation; the choice of these parameters depends on the
     * actual AMG version (subroutine crsng)
     */
    ecg1 = 0.0;
    ecg2 = 0.25;
    ewt2 = 0.35;

    /**
     * nwt - integer parameter affecting the creation of coarser grids and/or the definition
     * of the interpolation; the choice of this parameter depends on the actual AMG version
     * (subroutine crsng)
     */
    nwt = 2;

    /**
     * parameter controlling coarse-grid operator truncation
     * 0 - pairs of zeroes are removed from coarse grid operators
     * 1 - no coarse grid operator truncation
     */
    ntr = 0;
}

