RETISTRUCT_VERSION=$(shell grep Version retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz

dist: package
	rm -f retistruct.zip
	zip -r retistruct.zip $(PACKAGE) matlab  doc/retistruct-manual.pdf install.R -x '*/*.svn/*'

package:
	R CMD roxygen retistruct retistruct
	R CMD build retistruct

install:package
	R CMD INSTALL $(PACKAGE)