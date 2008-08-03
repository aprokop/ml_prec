AM_CPPFLAGS = -I$(top_srcdir) 

bin_PROGRAMS = main
main_SOURCES = main.cc main_helper.cc

main_LDADD = \
	$(top_builddir)/modules/libmodules.a 

all: prereq main 
	
prereq:
	$(MAKE) -C $(top_srcdir)/modules
