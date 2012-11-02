RETISTRUCT_VERSION=$(shell grep Version pkg/retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
RETISTRUCT_SVN_REVISION=$(shell svn info -R | grep "Revision:" | perl -p -e 's/Revision: //;' | sort -n -r | head -1)
RETISTRUCT_SVN_REVISION1=$(shell echo $(RETISTRUCT_SVN_REVISION) + 1 | bc) 
ifeq ("$(RETISTRUCT_SVN_REVISION)", "")
	RETISTRUCT_SVN_REVISION=1000
	RETISTRUCT_SVN_REVISION1=1000
endif
RETISTRUCT_PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz
RETISTRUCTGUI_PACKAGE=retistructgui_$(RETISTRUCT_VERSION).tar.gz
RETISTRUCTDEMOS_PACKAGE=retistructdemos_$(RETISTRUCT_VERSION).tar.gz

dist: retistruct retistructgui check user-guide
	rm -f retistruct_$(RETISTRUCT_VERSION).zip
	cp retistruct.Rcheck/retistruct-manual.pdf doc
	zip -r retistruct_$(RETISTRUCT_VERSION).zip $(RETISTRUCT_PACKAGE) $(RETISTRUCTGUI_PACKAGE) matlab  doc/retistruct-user-guide.pdf doc/retistruct-manual.pdf install.R install-gui.R README.txt  -x '*/*.svn/*' -x 'matlab/images/*'

roxygen:
	rm -f pkg/retistruct/man/*
	echo "library(roxygen2) ; roxygenize(\"pkg/retistruct\")" |	R --no-restore --slave

fix-revision:
	perl -p -i -e "s/^retistruct\.global\.revision.*/retistruct.global.revision <- $(RETISTRUCT_SVN_REVISION1)/;" pkg/retistruct/R/revision.R

retistruct: fix-revision roxygen 
	R CMD build pkg/retistruct

retistructgui: retistruct retistructdemos
	R CMD build pkg/retistructgui

retistructdemos: retistruct
	R CMD build pkg/retistructdemos

install: install-retistruct  install-retistructdemos install-retistructgui

install-retistruct: retistruct
	R CMD INSTALL --latex $(RETISTRUCT_PACKAGE) 

install-retistructgui: retistructgui
	R CMD INSTALL --latex $(RETISTRUCTGUI_PACKAGE) 

install-retistructdemos: retistructdemos
	R CMD INSTALL --latex $(RETISTRUCTDEMOS_PACKAGE) 

doc: retistruct
	rm -f retistruct.pdf
	R CMD Rd2dvi --pdf retistruct

user-guide:
	cd trunk/doc &&	pdflatex retistruct-user-guide.tex && cp retistruct-user-guide.pdf ../../www

check:
	R CMD check $(RETISTRUCT_PACKAGE)

revision:
	@echo $(RETISTRUCT_SVN_REVISION)
	@echo $(RETISTRUCT_SVN_REVISION1)

clean:
	rm -f pkg/retistruct/R/*~
	rm -f trunk/matlab/*~