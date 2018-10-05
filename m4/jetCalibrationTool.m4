#AC_SEARCH_JETCALIBRATIONTOOL(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_JETCALIBRATIONTOOL],[

  if test x$with_JetCalibrationTool != x && test x$with_JetCalibrationTool != xyes; then
    AC_MSG_NOTICE([Adding $with_JetCalibrationTool to search path for JetCalibrationTool])
    
    if test -d $with_JetCalibrationTool/ApplyJetCalibration && test -d $with_JetCalibrationTool/Standalone ; then
      found_JetCalibrationTool=yes
      JetCalibrationTool_include=$with_JetCalibrationTool/ApplyJetCalibration
      JetCalibrationTool_lib=$with_JetCalibrationTool/Standalone
    else
      found_JetCalibrationTool=no
    fi
  fi

  if test "x$found_JetCalibrationTool" = "xyes" ; then
    if test -f $JetCalibrationTool_include/ApplyJetCalibration.h && test -f $JetCalibrationTool_lib/libApplyJetCalibration.$LIB_SUFFIX ; then
      JETCALIBRATIONTOOL_LDFLAGS="-L$JetCalibrationTool_lib -lApplyJetCalibration"
      JETCALIBRATIONTOOL_CPPFLAGS="-I$JetCalibrationTool_include"
    else
      found_JetCalibrationTool = no
    fi
  fi

  if test x$found_JetCalibrationTool != xyes ; then
    for ac_jet_path_tmp in /usr /usr/local /opt /opt/local ; do
      if test -d $ac_jet_path_tmp/include/ApplyJetCalibration && test -d $ac_jet_path_tmp/lib; then
        found_JetCalibrationTool = yes
        JetCalibrationTool_include = $ac_jet_path_tmp/include/ApplyJetCalibration
        JetCalibrationTool_lib = $ac_jet_path_tmp/lib
        if test -f $JetCalibrationTool_include/ApplyJetCalibration.h && test -f $JetCalibrationTool_lib/libApplyJetCalibration.$LIB_SUFFIX ; then
          JETCALIBRATIONTOOL_LDFLAGS="-L$JetCalibrationTool_lib -lApplyJetCalibration"
          JETCALIBRATIONTOOL_CPPFLAGS="-I$JetCalibrationTool_include"
          break;
        else
          found_JetCalibrationTool = no
        fi
      fi
    done
  fi


  if test x$found_JetCalibrationTool = xyes ; then
  
    export JETCALIBRATIONTOOL_LDFLAGS
    export JETCALIBRATIONTOOL_CPPFLAGS
    AC_SUBST([JETCALIBRATIONTOOL_LDFLAGS])
    AC_SUBST([JETCALIBRATIONTOOL_CPPFLAGS])
    AC_MSG_NOTICE([Found JetCalibrationTool package])
    AC_MSG_NOTICE([JETCALIBRATIONTOOL_LDFLAGS = $JETCALIBRATIONTOOL_LDFLAGS])
    AC_MSG_NOTICE([JETCALIBRATIONTOOL_CPPFLAGS = $JETCALIBRATIONTOOL_CPPFLAGS])
    $1
  else
    AC_MSG_NOTICE([Could not find JetCalibrationTool package])
    $2
  fi

])