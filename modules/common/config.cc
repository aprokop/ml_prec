#include "common.h"
#include "include/exception.h"

std::ostream& operator<<(std::ostream& os, const Config& cfg) {
    os << "Level iterations : ";
    for (uint i = 0; i < cfg.niters.size(); i++)
	os << cfg.niters[i] << " ";
    os << std::endl;
    os << "Level sigmas     : ";
    for (uint i = 0; i < cfg.sigmas.size(); i++)
	os << cfg.sigmas[i] << " ";
    os << std::endl;

    /* Mesh parameters */
    if (cfg.matrix.empty()) {
	os << "c                : " << cfg.c << std::endl;
	os << "Geometry         : " << cfg.nx << " x " << cfg.ny << " x " << cfg.nz << std::endl;
	os << "Unsymmetric      : " << (cfg.unsym_matrix ? "true" : "false") << std::endl;
	if (cfg.unsym_matrix)
	    os << "Unsymmetrix shift: " << cfg.unsym_shift << std::endl;
    } else {
	os << "Matrix file      : " << cfg.matrix << std::endl;
    }
    os << "Dump data        : " << (cfg.dump_data ? "true" : "false") << std::endl;

    if (cfg.analysis != ANAL_NONE) {
	os << "Analysis         : ";
	switch(cfg.analysis) {
	    case ANAL_HISTOGRAMM: os << "histogramm"; break;
	    case ANAL_QDROPPED  : os << "qdropped"; break;
	    case ANAL_Q_REM_FIXED_ROW: os << "q_rem_fixed_row"; break;
	    default 		: THROW_EXCEPTION("Unknown AnalType");
	}
	os << std::endl;

	return os;
    }

    /* Run parameters */
    os << "Number of tests  : " << cfg.ntests << std::endl;
    os << "Use tails        : " << (cfg.use_tails ? "true" : "false") << std::endl;
    os << "Optimize storage : " << (cfg.optimize_storage ? "true" : "false") << std::endl;

    os << "Solver           : ";
    switch(cfg.solver) {
	case PCG_SOLVER    : os << "pcg"; break;
	case CHEB_SOLVER   : os << "cheb"; break;
	case SIMPLE_SOLVER : os << "simple"; break;
	case DIRECT_SOLVER : os << "direct"; break;
	default		: THROW_EXCEPTION("Unknown SolverType");
    }
    os << std::endl;

    os << "Preconditioner   : ";
    switch(cfg.prec) {
	case UH_CHEB_PREC     : os << "uh_cheb"; break;
	case AMG_PREC         : os << "amg"; break;
	case DIAG_PREC        : os << "diag"; break;
	case GS_PREC          : os << "gs"; break;
	case SYM_SPLIT_PREC   : os << "sym_split"; break;
	case MULTI_SPLIT_PREC : os << "multi_split"; break;
	default           : THROW_EXCEPTION("Unknown PrecType");
    }
    os << std::endl;

    return os;
}
