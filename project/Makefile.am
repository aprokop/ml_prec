AM_CPPFLAGS = -I$(top_srcdir)

if HYPRE
    HYP=hypre
else
    HYP=
endif

bin_PROGRAMS = ml_prec $(HYP) convert
ml_prec_SOURCES = main.cc main_helper.cc analysis.cc transform.cc main.h config.h
hypre_SOURCES = hypre.cc

ml_prec_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lrt

hypre_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lHYPRE -lrt

convert_SOURCES = convert.cc
convert_LDADD = $(top_builddir)/modules/libmodules.a

all: prereq ml_prec $(HYP)

prereq:
	$(MAKE) -C $(top_srcdir)/modules
