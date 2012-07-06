RETISTRUCT_VERSION=$(shell grep Version retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
RETISTRUCT_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
RETISTRUCT_SVN_REVISION1=$(shell echo $(RETISTRUCT_SVN_REVISION) + 1 | bc) 
RETISTRUCT_PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz
RETISTRUCTGUI_PACKAGE=retistructgui_$(RETISTRUCT_VERSION).tar.gz

dist: retistruct retistructgui
	rm -f retistruct_$(RETISTRUCT_VERSION).zip
	zip -r retistruct_$(RETISTRUCT_VERSION).zip $(RETISTRUCT_PACKAGE) $(RETISTRUCTGUI_PACKAGE) matlab  doc/retistruct-manual.pdf install.R install-gui.R README.txt -x '*/*.svn/*' -x 'matlab/images/*'

roxygen:
	rm -f retistruct/man/*
	echo "library(roxygen2) ; roxygenize(\"retistruct\")" |	R --no-restore --slave

fix-revision:
	perl -p -i -e "s/^retistruct\.global\.revision.*/retistruct.global.revision <- $(RETISTRUCT_SVN_REVISION1)/;" retistruct/R/revision.R

retistruct: fix-revision roxygen 
	R CMD build retistruct

retistructgui: retistruct
	R CMD build retistructgui

install: install-retistruct install-retistructgui

install-retistruct: retistruct
	R CMD INSTALL --latex $(RETISTRUCT_PACKAGE) 

install-retistructgui: retistructgui
	R CMD INSTALL --latex $(RETISTRUCTGUI_PACKAGE) 

doc: retistruct
	rm -f retistruct.pdf
	R CMD Rd2dvi --pdf retistruct

check:
	R CMD check $(RETISTRUCT_PACKAGE)

revision:
	@echo $(RETISTRUCT_SVN_REVISION)
	@echo $(RETISTRUCT_SVN_REVISION1)

clean:
	rm -f retistruct/R/*~
	rm -f matlab/*~