AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)

noinst_LIBRARIES = libmodules.a
libmodules_a_SOURCES = \
		       config_file/config_file.cc \
		       config_file/config_file.h \
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
		       prec/amg/amg_prec.cc \
		       prec/amg/amg_prec.h \
		       prec/cheb/cheb_constr.cc \
		       prec/cheb/cheb_mark.cc \
		       prec/cheb/cheb_misc.cc \
		       prec/cheb/cheb_solve.cc \
		       prec/cheb/cheb_misc.h \
		       prec/cheb/cheb_prec.h \
		       prec/composite/composite_prec.cc \
		       prec/composite/composite_prec.h \
		       prec/diag/diag_prec.cc \
		       prec/diag/diag_prec.h \
		       prec/gs/gs_prec.cc \
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
		       prec/sym/sym_prec.cc \
		       prec/sym/sym_prec.h \
		       prec/prec_base.h \
		       prec/prec.h \
		       solvers/cheb.cc \
		       solvers/pcg.cc \
		       solvers/simple_solver.cc \
		       solvers/gmres.cc \
		       solvers/amg1r6.f \
		       solvers/solvers.cc \
		       solvers/solvers.h
if UMFPACK
libmodules_a_SOURCES += \
		       prec/bgs/bgs_prec.cc \
		       prec/bgs/bgs_prec.h \
		       prec/rbgs/rbgs_prec.cc \
		       prec/rbgs/rbgs_prec.h \
		       solvers/umfpack.cc
endif
