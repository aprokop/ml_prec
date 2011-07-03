AM_CPPFLAGS = -I$(top_srcdir)

noinst_LIBRARIES = libmodules.a
libmodules_a_SOURCES = \
		       matrix/csrmatrix.cc \
		       matrix/dmatrix.cc \
		       matrix/mapmatrix.cc \
		       matrix/skyline.cc \
		       matrix/matrix.h \
		       vector/blas.cc \
		       vector/vector.cc \
		       vector/vector.h \
		       mesh/spe.cc \
		       mesh/mesh.h \
		       prec/amg/amg.cc \
		       prec/amg/amg_prec.h \
		       prec/cheb/cheb_constr.cc \
		       prec/cheb/cheb_mark.cc \
		       prec/cheb/cheb_misc.cc \
		       prec/cheb/cheb_solve.cc \
		       prec/cheb/cheb_misc.h \
		       prec/cheb/cheb_prec.h \
		       prec/diag/diag.cc \
		       prec/diag/diag_prec.h \
		       prec/gs/gs.cc \
		       prec/gs/gs_prec.h \
		       prec/id/id_prec.h \
		       prec/misc/misc.cc \
		       prec/misc/sparse_lu.cc \
		       prec/misc/misc.h \
		       prec/multi_split/multi_split_constr.cc \
		       prec/multi_split/multi_split_mark.cc \
		       prec/multi_split/multi_split_misc.cc \
		       prec/multi_split/multi_split_solve.cc \
		       prec/multi_split/multi_split_misc.h \
		       prec/multi_split/multi_split_prec.h \
		       prec/sym/sym.cc \
		       prec/sym/sym_prec.h \
		       prec/prec_base.h \
		       prec/prec.h \
		       solvers/cheb.cc \
		       solvers/pcg.cc \
		       solvers/simple_solver.cc \
		       solvers/amg1r6.f \
		       solvers/solvers.cc \
		       solvers/solvers.h
if UMFPACK
libmodules_a_SOURCES += \
		       prec/bgs/bgs.cc \
		       prec/bgs/bgs_prec.h \
		       prec/rbgs/rbgs.cc \
		       prec/rbgs/rbgs_prec.h \
		       solvers/umfpack.cc
endif
