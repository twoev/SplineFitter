#AC_SEARCH_GOODRUNSLISTS(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_GOODRUNSLISTS],[
  
  if test x$with_GoodRunsLists != x && test x$with_GoodRunsLists != xyes ; then
    AC_MSG_NOTICE([Adding $with_GoodRunsLists to search path for GoodRunsLists])
    if test -d $with_GoodRunsLists/include/GoodRunsLists && test -d $with_GoodRunsLists/lib ; then
      found_GoodRunsLists=yes
      GoodRunsLists_include=$with_GoodRunsLists/include
      GoodRunsLists_lib=$with_GoodRunsLists/lib
    else
      found_GoodRunsLists = no
    fi
  fi
  
  if test "x$found_GoodRunsLists" = "xyes" ; then
    if test -f $GoodRunsLists_include/GoodRunsLists/TGoodRunsListReader.h && test -f $GoodRunsLists_include/GoodRunsLists/TGoodRunsList.h ; then
      GOODRUNSLISTS_LDFLAGS="-L$GoodRunsLists_lib -lGoodRunsLists"
      GOODRUNSLISTS_CPPFLAGS="-I$GoodRunsLists_include"
    else
    found_GoodRunsList = no
    fi
  fi
  
  if test x$found_GoodRunsLists != xyes; then
    for ac_grl_path_tmp in /usr /usr/local /opt /opt/local ; do
      if test -d $ac_grl_path_tmp/include/GoodRunsLists && test -d $ac_grl_path_tmp/lib; then
        found_GoodRunsLists = yes
        GoodRunsLists_include = $ac_grl_path_tmp/include/GoodRunsLists
        GoodRunsLists_lib = $ac_grl_path_tmp/lib
        if test -f $GoodRunsLists_include/TGoodRunsListReader.h && test -f $GoodRunsLists_include/TGoodRunsList.h ; then
          GOODRUNSLISTS_LDFLAGS = "-L$GoodRunsLists_lib -lGoodRunsLists"
          GOODRUNSLISTS_CPPFLAGS = "-I$GoodRunsLists_include"
          break;
        else
          found_GoodRunsLists = no
        fi  
      fi
    done
  fi
  
  if test x$found_GoodRunsLists = xyes ; then
  
    export GOODRUNSLISTS_LDFLAGS
    export GOODRUNSLISTS_CPPFLAGS
    AC_SUBST([GOODRUNSLISTS_LDFLAGS])
    AC_SUBST([GOODRUNSLISTS_CPPFLAGS])
    AC_MSG_NOTICE([Found GoodRunsLists package])
    AC_MSG_NOTICE([GOODRUNSLISTS_LDFLAGS = $GOODRUNSLISTS_LDFLAGS])
    AC_MSG_NOTICE([GOODRUNSLISTS_CPPFLAGS = $GOODRUNSLISTS_CPPFLAGS])
    $1
  else
    AC_MSG_NOTICE([Could not find GoodRunsLists package])
    $2
  fi
  
])

