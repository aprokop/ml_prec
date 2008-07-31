AM_CPPFLAGS = -I$(top_srcdir) 

bin_PROGRAMS = _main
_main_SOURCES = main.cc 

_main_LDADD = \
	$(top_builddir)/modules/libmodules.a 

all: prereq _main 
	
prereq:
	$(MAKE) -C $(top_srcdir)/modules
