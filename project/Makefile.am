AM_CPPFLAGS = -I$(top_srcdir) 

bin_PROGRAMS = spe_prec
spe_prec_SOURCES = main.cc main_helper.cc

spe_prec_LDADD = \
	$(top_builddir)/modules/libmodules.a \
	-lrt

all: prereq spe_prec 
	
prereq:
	$(MAKE) -C $(top_srcdir)/modules
