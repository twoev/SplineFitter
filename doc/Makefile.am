if DX_COND_doc

.PHONY: doxygen-run doxygen-doc clean

.INTERMEDIATE: doxygen-run

doxygen-run: @DX_DOCDIR@/@PACKAGE@.tag

doxygen-doc: doxygen-run 

doc: doxygen-doc

@DX_DOCDIR@/@PACKAGE@.tag: $(DX_CONFIG) $(pkginclude_HEADERS)
	rm -rf @DX_DOCDIR@
	$(DX_ENV) $(DX_DOXYGEN) $(srcdir)/$(DX_CONFIG)

clean: 
	rm -rf html latex

endif DX_COND_doc