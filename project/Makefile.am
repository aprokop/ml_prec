AM_CPPFLAGS = -I$(top_srcdir)

if HYPRE
    HYP=hypre
else
    HYP=
endif

bin_PROGRAMS = spe_prec $(HYP) convert
spe_prec_SOURCES = main.cc main_helper.cc analysis.cc
hypre_SOURCES = hypre.cc

spe_prec_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lumfpack \
	-lrt

hypre_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lHYPRE -lrt

convert_SOURCES = convert.cc
convert_LDADD = $(top_builddir)/modules/libmodules.a

all: prereq spe_prec $(HYP)

prereq:
	$(MAKE) -C $(top_srcdir)/modules
