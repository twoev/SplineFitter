#AC_SEARCH_GSL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_GSL],[
GSLCONFIGPATH=$PATH
if test x$with_gsl != x; then
  AC_MSG_NOTICE([Adding $with_gsl to gsl-config search path])
  GSLCONFIGPATH=$with_gsl/bin:$GSLCONFIGPATH
fi

AC_MSG_NOTICE([gsl-config search path is $GSLCONFIGPATH])

AC_PATH_PROG(GSLCONFIG, [gsl-config], [no], [$GSLCONFIGPATH])

if test x$GSLCONFIG != xno; then
  AC_MSG_NOTICE([Found gsl-config])
  GSLCPPFLAGS="`$GSLCONFIG --cflags`"
  GSLLIBS="`$GSLCONFIG --libs`"

  AC_SUBST([GSLCPPFLAGS])
  AC_SUBST([GSLLIBS])

  AC_MSG_NOTICE([GSLCPPFLAGS=$GSLCPPFLAGS])
  AC_MSG_NOTICE([GSLLIBS=$GSLLIBS])

  $1
else
  AC_MSG_NOTICE([GNU Scientific library not found!!])
  $2
fi

])