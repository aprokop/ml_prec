# Prelude
AC_INIT([ml_prec], [0.9.6], [andrey.prok@gmail.com])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

AM_PROG_AR([ar])

AC_LANG([C++])

# Check for programs
AC_PROG_CXX([g++])
AC_PROG_F77
AC_PROG_LIBTOOL
AC_PROG_LN_S

AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN

# Check for header files
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Check for libraries

# Optimization
FFLAGS="-O3"
CFLAGS="-O3"
CXXFLAGS="-W -Wall -Wno-unused"
# CXXFLAGS="-W -Wall -Wno-unused -pedantic -ansi -pthread"

# DEBUG
AC_MSG_CHECKING([wish for debug])
AC_ARG_WITH(debug,
            [AC_HELP_STRING([--with-debug@<:@=full|assert|none@:>@],
                            [enable debugging and parameters check @<:@default=none@:>@])],
            __debug=$withval,
            __debug="none")
AC_MSG_RESULT($__debug)
test x$__debug = xfull   && CXXFLAGS="-g -O0 $CXXFLAGS"
test x$__debug = xassert && CXXFLAGS="-O2 $CXXFLAGS"
test x$__debug = xnone   && CXXFLAGS="-DNO_ASSERT -O3 $CXXFLAGS"

LIBS="-lm $LIBS"

# Boost
AC_CHECK_HEADER([boost/version.hpp],
                [AC_DEFINE([HAVE_BOOST], [1], [Define to 1 if you have the Boost library])],
                [])

# BLAS
AC_MSG_CHECKING([wish for blas])
AC_ARG_ENABLE(blas,
              [AC_HELP_STRING([--enable-blas],
                              [use BLAS @<:@default=yes@:>@])],
              enable_blas=$enableval,
              enable_blas="yes")
AC_MSG_RESULT($enable_blas)
if test x$enable_blas = xyes; then
    AX_BLAS
fi
LIBS="$BLAS_LIBS $LIBS $FLIBS"

# UMFPACK
AC_CHECK_HEADERS([umfpack.h suitesparse/umfpack.h], [enable_umfpack="yes"], [enable_umfpack="no"])
if test x$enable_umfpack = xyes; then
    AC_DEFINE([HAVE_UMFPACK], [1], [Define to 1 if you have UMFPACK library])
    LIBS="-lumfpack $LIBS"
fi
AM_CONDITIONAL([UMFPACK], [test x$enable_umfpack = xyes])

# Logger
AC_MSG_CHECKING([wish for logger])
AC_ARG_ENABLE(logger,
              [AC_HELP_STRING([--enable-logger],
                              [use Apache log4cxx @<:@default=yes@:>@])],
              enable_logger=$enableval,
              enable_logger="yes")
AC_MSG_RESULT($enable_logger)
if test x$enable_logger = xyes; then
    AC_CHECK_HEADER([log4cxx/propertyconfigurator.h],
                    [LIBS="-llog4cxx $LIBS"],
                    [AC_DEFINE([NO_LOGGER], [1], [Disable log4cxx])])
else
    AC_DEFINE([NO_LOGGER], [1], [Disable logger])
fi

# ViTE
AC_MSG_CHECKING([wish for ViTE])
AC_ARG_ENABLE(vite,
              [AC_HELP_STRING([--enable-vite],
                              [use ViTE tracer @<:@default=yes@:>@])],
              enable_vite=$enableval,
              enable_vite="yes")
AC_MSG_RESULT($enable_vite)
if test x$enable_vite = xyes; then
    AC_DEFINE([USE_VITE], [1], [Output tracing info in Paje format])
fi

# HYPRE
AC_MSG_CHECKING([wish for HYPRE])
AC_ARG_ENABLE(hypre,
              [AC_HELP_STRING([--enable-hypre],
                      [enable HYPRE framework @<:@default=no@:>@])],
              enable_hypre=$withval,
              enable_hypre="no")
AC_MSG_RESULT($enable_hypre)
AM_CONDITIONAL([HYPRE], [test x$enable_hypre = xyes])
if test x$enable_hypre = xyes; then
    AX_MPI
    LIBS="$LIBS $MPILIBS"
fi

# gprof
AC_MSG_CHECKING([wish for gprof])
AC_ARG_ENABLE(profiling,
              [AC_HELP_STRING([--enable-profiling],
                              [enable gprof @<:@default=no@:>@])],
              enable_prof=$enableval,
              enable_prof="no")
AC_MSG_RESULT($enable_prof)
if test x$enable_prof = xyes; then
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
