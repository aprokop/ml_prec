AM_CPPFLAGS = -I$(top_srcdir)

if HYPRE
    HYP=hypre
else
    HYP=
endif

bin_PROGRAMS = spe_prec $(HYP)
spe_prec_SOURCES = main.cc main_helper.cc
hypre_SOURCES = hypre.cc

spe_prec_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lrt

hypre_LDADD = -lHYPRE -lrt

all: prereq spe_prec $(HYP)

prereq:
	$(MAKE) -C $(top_srcdir)/modules
