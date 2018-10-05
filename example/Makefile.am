SUBDIRS = . tclap

bin_PROGRAMS = FitBSpline

FitBSpline_SOURCES = main.cxx

FitBSpline_CPPFLAGS  = $(AM_CPPFLAGS) -I$(top_builddir) $(ROOTCPPFLAGS)
FitBSpline_LDFLAGS = $(AM_LDFLAGS) $(ROOTLDFLAGS)
FitBSpline_LDADD = $(AM_LDADD) $(top_builddir)/src/.libs/libSpline.la  $(ROOTLIBS)

