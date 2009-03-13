# Prelude
AC_INIT([spe_prec], [0.1.0])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_LANG([C++])

# Check for programs
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77([gfortran])
AC_PROG_LIBTOOL
AC_PROG_LN_S

# Check for libraries

# FFLAGS
FFLAGS="-O3"
# CXXFLAGS
CFLAGS="-O3"
CXXFLAGS="-W -Wall -Wno-unused -pedantic -ansi -pthread"
AC_ARG_WITH(debug, AC_HELP_STRING([--with-debug@<:@=full|assert|none@:>@], [enable debugging and parameters check @<:@default=assert@:>@]),__debug="$withval", __debug="assert")
test "x$__debug" = "xfull"   && CXXFLAGS="-g -O0 $CXXFLAGS"
test "x$__debug" = "xassert" && CXXFLAGS="-O3 $CXXFLAGS"
test "x$__debug" = "xnone"   && CXXFLAGS="-DNO_ASSERT -O3 $CXXFLAGS"

# Check for header files
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

LIBS="-lm $LIBS"

# Boost
AC_CHECK_HEADER([boost/version.hpp], [AC_DEFINE([HAVE_BOOST], [1], [Define to 1 if you have the Boost library])], [])

# Logger
AC_CHECK_HEADER( [log4cxx/propertyconfigurator.h], [LIBS="-llog4cxx $LIBS"], [AC_DEFINE(NO_LOGGER, [1], [Disable logger])] )

# BLAS
AC_ARG_WITH(blas, AC_HELP_STRING([--with-blas@<:@=yes|no@:>@], [use BLAS @<:@default=yes@:>@]), __blas="$withval", __blas="yes")
test "x$__blas" = "xyes" && ACX_BLAS
LIBS="$BLAS_LIBS $LIBS"

# gprof
AC_ARG_WITH(profiling, AC_HELP_STRING([--with-profiling@<:@=yes|no@:>@], [enable gprof @<:@default=no@:>@]), __prof="$withval", __prof="no")
if test "x$__prof" = "xyes"; then
    CXXFLAGS="-pg $CXXFLAGS"
    LDFLAGS="-pg $LDFLAGS"
fi

# Output files
AC_CONFIG_HEADERS(config/config.h)
AC_CONFIG_FILES(
	Makefile
	modules/Makefile 
	project/Makefile
	tools/Makefile
	tools/check_memory/Makefile
)

AC_OUTPUT
