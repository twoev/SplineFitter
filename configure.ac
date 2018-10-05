AC_PREREQ(2.59)

AC_INIT( [Spline], [0.1.0], [jmonk@cern.ch], [Spline])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_SRCDIR([Spline/BSpline.hh])
AC_CONFIG_HEADERS([config/config.h])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])
AC_SUBST(LT_OBJDIR)
AC_LANG(C++)

AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S

DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([Spline],, [doc])

## This macro checks for boost in autoconf 2.59
AX_BOOST_BASE([1.33.0], [AC_MSG_NOTICE([Found Boost])], [AC_MSG_ERROR([Boost library cannot be found!])])

# lib suffix differs on OS X
AC_DEFINE_DYNAMIC_LIBS

AC_ARG_ENABLE([root],
  [AC_HELP_STRING(--disable-root, [Do not link against ROOT - the example will not be built (default: enabled)])],
  [], [enable_root=yes])

if test x$enable_root = xyes; then
  AC_SEARCH_ROOT([AC_MSG_NOTICE([Found ROOT])], [AC_MSG_ERROR([ROOT is required and is not found])])
##These are for tclap
  AC_CXX_HAVE_SSTREAM
  AC_CXX_HAVE_STRSTREAM
else
  AC_MSG_NOTICE([*** Not building with ROOT - the example program will not be enabled!! ***])
fi

AM_CONDITIONAL(ENABLE_ROOT, [test x$enable_root = xyes])

## GNU Scientific Library

AC_ARG_WITH([gsl],
  [AC_HELP_STRING([--with-gsl], [Give the GSL install location])],
  )
  
AC_SEARCH_GSL([], [AC_MSG_ERROR([GSL requested but not found!!])])
  
AM_CXXFLAGS="$AM_CXXFLAGS -Wall -O2"
AM_CPPFLAGS="$AM_CPPFLAGS -Wall -O2 -I\$(top_builddir) $BOOST_CPPFLAGS"
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_CPPFLAGS)

AC_CONFIG_FILES([ Makefile 
                  src/Makefile
                  Spline/Makefile
                  example/Makefile
                  example/tclap/Makefile
                  doc/Makefile ])
                  
AC_OUTPUT
