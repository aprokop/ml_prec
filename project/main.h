#ifndef __MAIN_H__
#define __MAIN_H__

#include "include/define.h"
#include "modules/mesh/mesh.h"
#include "modules/matrix/matrix.h"
#include "modules/common/common.h"
#include <iostream>
#include <vector>

int set_params(int argc, char * argv[], Config& config);
double avg_time(const std::vector<double>& times);

void construct_matrix(const Config& cfg, const SPEMesh& mesh, SkylineMatrix& A);
void construct_vector(const Config& cfg, Vector& b);
void dump_data(const SkylineMatrix& A, const Vector& b);

void examine(const SkylineMatrix& A, AnalType analysis);
void analysis(const SkylineMatrix& A);

struct GlobalStats {
    double t_resid;	    /* Residual calculation time */
    double t_prec;	    /* Preconditioner inversion time */
    double t_const;	    /* Preconditioner construction time */
    double t_sol;	    /* System solution time */
    double t_total;	    /* = t_const + t_sol */
    uint   niter;	    /* Number of iterations */

    void dump(const std::string& dir) const;
};

#endif // __MAIN_H__
