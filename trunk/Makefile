RETISTRUCT_VERSION=$(shell grep Version retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
RETISTRUCT_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
RETISTRUCT_SVN_REVISION1=$(shell echo $(RETISTRUCT_SVN_REVISION) + 1 | bc) 
PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz

dist: package
	rm -f retistruct.zip
	zip -r retistruct.zip $(PACKAGE) matlab  doc/retistruct-manual.pdf install.R -x '*/*.svn/*'

roxygen:
	rm -f retistruct/man/*
	R CMD roxygen -d retistruct

fix-revision:
	perl -p -i -e "s/^retistruct\.global\.revision.*/retistruct.global.revision <- $(RETISTRUCT_SVN_REVISION1)/;" retistruct/R/revision.R

package: fix-revision roxygen 
	R CMD build retistruct

install: package
	R CMD INSTALL --latex $(PACKAGE) 

doc: package
	rm -f retistruct.pdf
	R CMD Rd2dvi --pdf retistruct

check:
	R CMD check $(PACKAGE)

revision:
	@echo $(RETISTRUCT_SVN_REVISION)
	@echo $(RETISTRUCT_SVN_REVISION1)

clean:
	rm -f retistruct/R/*~