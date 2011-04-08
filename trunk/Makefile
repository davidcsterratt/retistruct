RETISTRUCT_VERSION=$(shell grep Version retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
RETISTRUCT_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz



dist: package
	rm -f retistruct.zip
	zip -r retistruct.zip $(PACKAGE) matlab  doc/retistruct-manual.pdf install.R -x '*/*.svn/*'

package:
	perl -p -i -e "s/^retistruct\.global\.revision.*/retistruct.global.revision <- $(RETISTRUCT_SVN_REVISION)/;" retistruct/R/retistruct.R
	R CMD roxygen retistruct retistruct
	R CMD build retistruct

install:package
	R CMD INSTALL $(PACKAGE)

revision:
	echo $(RETISTRUCT_SVN_REVISION)
