#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

#include "include/uvector.h"
#include "include/time.h"

#include <cmath>
#include <cstdlib>
#include <mpi.h>

int main (int argc, char *argv[]) {
    /* Initialize MPI */
    MPI_Init(&argc, &argv);

    /* Read command-line parameters */
    int solver_id = 0;
    while (1) {
	int option_index = 0;
	int ch = getopt(argc, argv, "hs:");

	if (ch == -1)
	    break;

	if (ch == 0)
	    continue;

	switch (ch) {
	    case 'h':
		std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "  -s {0|1|8|50|61}        Solver ID (AMG | AMG-PCG | ParaSails-PCG | PCG | AMG-FlexGMRES" << std::endl;
		return 0;
	    case 's':
		solver_id = atoi(optarg);
		break;
	    case '?':
	    default :
		abort();
	}
    }

    /* Read the matrix from file */
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;

    HYPRE_IJMatrixRead("matrix_hypre.dat", MPI_COMM_WORLD, HYPRE_PARCSR, &A);
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

    int ilower, iupper, jlower, jupper;
    HYPRE_IJMatrixGetLocalRange(A, &ilower, &iupper, &jlower, &jupper);

    /* Read the vector from file */
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;

    HYPRE_IJVectorRead("vector_hypre.dat", MPI_COMM_WORLD, HYPRE_PARCSR, &b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);

    /* Set space for x */
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    /* Set x to random */
    int local_size = iupper - ilower + 1;
    uvector<double> xval(local_size);
    uvector<int> rows(local_size);
    srandom(3);
    for (int i = 0; i < local_size; i++) {
	xval[i] = 20.*(random() - 0.5*RAND_MAX)/RAND_MAX + 100;
	rows[i] = i;
    }
    HYPRE_IJVectorSetValues(x, local_size, &rows[0], &xval[0]);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    HYPRE_Solver solver, precond;

    /* Choose a solver and solve the system */
    double ctime = 0, stime = 0;
    /* AMG */
    if (solver_id == 0)
    {
	int num_iterations;
	double final_res_norm;

	/* Create solver */
	HYPRE_BoomerAMGCreate(&solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
	HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
	HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
	HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
	HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
	HYPRE_BoomerAMGSetTol(solver, 1e-7);      /* conv. tolerance */

	/* Now setup and solve! */
	ctime = pclock();
	HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
	ctime = pclock() - ctime;

	stime = pclock();
	HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
	stime = pclock() - stime;

	/* Run info - needed logging turned on */
	HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
	HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destroy solver */
	HYPRE_BoomerAMGDestroy(solver);
    }
    /* PCG */
    else if (solver_id == 50)
    {
	int num_iterations;
	double final_res_norm;

	/* Create solver */
	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
	HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
	HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
	HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

	/* Now setup and solve! */
	HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

	/* Run info - needed logging turned on */
	HYPRE_PCGGetNumIterations(solver, &num_iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destroy solver */
	HYPRE_ParCSRPCGDestroy(solver);
    }
    /* PCG with AMG preconditioner */
    else if (solver_id == 1)
    {
	int num_iterations;
	double final_res_norm;

	/* Create solver */
	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
	HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
	HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

	/* Now set up the AMG preconditioner and specify any parameters */
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

	/* Set the PCG preconditioner */
	HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);

	/* Now setup and solve! */
	ctime = pclock();
	HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	ctime = pclock() - ctime;

	stime = pclock();
	HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);
	stime = pclock() - stime;

	/* Run info - needed logging turned on */
	HYPRE_PCGGetNumIterations(solver, &num_iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destroy solver and preconditioner */
	HYPRE_ParCSRPCGDestroy(solver);
	HYPRE_BoomerAMGDestroy(precond);

    } else if (solver_id == 8) {
	/* PCG with Parasails Preconditioner */
	int    num_iterations;
	double final_res_norm;

	int      sai_max_levels = 1;
	double   sai_threshold = 0.1;
	double   sai_filter = 0.05;
	int      sai_sym = 1;

	/* Create solver */
	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
	HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
	HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

	/* Now set up the ParaSails preconditioner and specify any parameters */
	HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
	HYPRE_ParaSailsSetFilter(precond, sai_filter);
	HYPRE_ParaSailsSetSym(precond, sai_sym);
	HYPRE_ParaSailsSetLogging(precond, 3);

	/* Set the PCG preconditioner */
	HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
			    (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);

	/* Now setup and solve! */
	HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);


	/* Run info - needed logging turned on */
	HYPRE_PCGGetNumIterations(solver, &num_iterations);
	HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destory solver and preconditioner */
	HYPRE_ParCSRPCGDestroy(solver);
	HYPRE_ParaSailsDestroy(precond);
    } else if (solver_id == 61) {
	/* Flexible GMRES with  AMG Preconditioner */
	int    num_iterations;
	double final_res_norm;
	int    restart = 30;
	int    modify = 1;


	/* Create solver */
	HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_FlexGMRESSetKDim(solver, restart);
	HYPRE_FlexGMRESSetMaxIter(solver, 1000); /* max iterations */
	HYPRE_FlexGMRESSetTol(solver, 1e-7); /* conv. tolerance */
	HYPRE_FlexGMRESSetPrintLevel(solver, 2); /* print solve info */
	HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */


	/* Now set up the AMG preconditioner and specify any parameters */
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

	/* Set the FlexGMRES preconditioner */
	HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
				  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);


	/* Now setup and solve! */
	HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
	HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

	/* Run info - needed logging turned on */
	HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
	HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

	printf("\n");
	printf("Iterations = %d\n", num_iterations);
	printf("Final Relative Residual Norm = %e\n", final_res_norm);
	printf("\n");

	/* Destory solver and preconditioner */
	HYPRE_ParCSRFlexGMRESDestroy(solver);
	HYPRE_BoomerAMGDestroy(precond);

    }
    else {
	std::cout << "Invalid solver id specified." << std::endl;
    }

    std::cout << "Construction time: " << ctime << std::endl;
    std::cout << "Solution time    : " << stime << std::endl;

    /* Print the solution */
#if 0
    HYPRE_IJVectorPrint(x, "ij.out.x");
#endif

    /* Clean up */
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);

    /* Finalize MPI*/
    MPI_Finalize();

    return(0);
}
