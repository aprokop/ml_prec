#include "cheb_prec.h"
#include "modules/prec/misc/misc.h"
#include "include/logger.h"

#include <map>

DEFINE_LOGGER("Prec");

void Prec::graph_planes(const std::string& filename, uint level, char plane, const SPEMesh& mesh) const {
    THROW_EXCEPTION("Not implemented");
#if 0
    // construct reverse map
    std::map<uint,uint> rev_map;
    if (level) {
	uvector<uint> gtr = levels[level-1].tr;
	uint n = gtr.size();
	for (int l = level-2; l >= 0; l--)
	    for (uint i = 0; i < n; i++)
		gtr[i] = levels[l].tr[gtr[i]];
	for (uint i = 0; i < n; i++)
	    rev_map[gtr[i]] = i;
	gtr.clear();

	::graph_planes(filename, levels[level].A, rev_map, plane, level, mesh);
    } else {
	uint n = level0_A.size();
	for (uint k = 0; k < n; k++)
	    rev_map[k] = k;

	::graph_planes(filename, level0_A, rev_map, plane, level, mesh);
    }
#endif
}

std::ostream& operator<<(std::ostream& os, const Prec& p) {
    os << std::endl;
    os << "nlevels = " << p.nlevels << std::endl;
    float oc = 0, gc = 0;
    for (uint level = 0; level < p.nlevels; level++) {
	gc += p.levels[level].N;
	oc += p.levels[level].nnz;
    }
    gc /= p.levels[0].N;
    oc /= p.levels[0].nnz;
    os << "gc     = " << gc << std::endl;   /* Grid complexity */
    os << "oc     = " << oc << std::endl;   /* Operator complexity */
    for (uint level = 0; level < p.nlevels; level++) {
	os << std::endl << "================== Level: " << level << " =======================" << std::endl;
	os << p.levels[level];
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const Prec::Level& li) {
    os << "N = " << li.N << ", nnz = " << li.nnz << std::endl;
    os << "M = " << li.M << ", Md  = " << li.Md << std::endl;
    os << "ncheb = " << li.ncheb << std::endl;
    os << "alpha = " << li.alpha << ", beta = " << li.beta << std::endl;
    os << "[lmin, lmax] = [" << li.lmin << "," << li.lmax << "], cond = " << li.lmax/li.lmin << std::endl;
    os << std::endl;
#if 0
    os << "A: " << li.A;
    os << "L: " << li.L;
    os << "U: " << li.U;
#endif
    return os;
}

void print_vite_header(std::ofstream& fos);
void Prec::dump_vite_trace() const {
    std::ofstream fos("vite.trace");
    print_vite_header(fos);

    fos << std::fixed << std::setprecision(8);

    fos << "7 0 T0 T TTP \"Thread 0\"" << std::endl;
    fos << (*foss).str();

    fos.close();

    // delete foss;
}

void print_vite_header(std::ofstream& fos) {
    fos << "\
%EventDef PajeDefineContainerType 1\n\
% Alias string\n\
% ContainerType string\n\
% Name string\n\
%EndEventDef\n\
%EventDef PajeDefineStateType 3\n\
% Alias string\n\
% ContainerType string\n\
% Name string\n\
%EndEventDef\n\
%EventDef PajeDefineEventType 4\n\
% Alias string\n\
% ContainerType string\n\
% Name string\n\
%EndEventDef\n\
%EventDef PajeDefineEntityValue 6\n\
% Alias string\n\
% EntityType string\n\
% Name string\n\
% Color color\n\
%EndEventDef\n\
%EventDef PajeCreateContainer 7\n\
% Time date\n\
% Alias string\n\
% Type string\n\
% Container string\n\
% Name string\n\
%EndEventDef\n\
%EventDef PajeDestroyContainer 8\n\
% Time date\n\
% Name string\n\
% Type string\n\
%EndEventDef\n\
%EventDef PajeSetState 10\n\
% Time date\n\
% Type string\n\
% Container string\n\
% Value string\n\
%EndEventDef\n\
%EventDef PajeNewEvent 11\n\
% Time date\n\
% Type string\n\
% Container string\n\
% Value int\n\
%EndEventDef\n\
1 P 0 Program\n\
1 T P Thread\n\
3 S T 'Thread State'\n\
4 E T 'Level'\n\
6   d	    S   'Default'	    '0.0 0.0 0.0'\n\
6   MV	    S   'MatrixVector'	    '1.0 0.0 0.0'\n\
6   Copy    S	'Copy'		    '0.5 0.5 0.5'\n\
6   Diag    S	'Diag'		    '0.0 1.0 0.0'\n\
6   LU	    S	'LU+nA'		    '0.0 0.5 0.5'\n\
6   M	    S	'Marking'	    '0.5 0.0 0.5'\n\
7 0 TTP P 0 'Program'\n";
}
