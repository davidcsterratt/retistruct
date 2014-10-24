RETISTRUCT_VERSION=$(shell grep Version pkg/retistruct/DESCRIPTION | perl -p -e "s/Version: //;")
RETISTRUCTDEMOS_VERSION=$(shell grep Version pkg/retistructdemos/DESCRIPTION | perl -p -e "s/Version: //;")
# RETISTRUCT_SVN_VERSION=$(shell git svn info | grep "Last Changed Rev:" | perl -p -e 's/Last Changed Rev: //;' | sort -n -r | head -1)
RETISTRUCT_GIT_COMMIT_HASH=$(shell git log --pretty="%h" HEAD^..HEAD)
RETISTRUCT_GIT_AUTHORDATE=$(shell git log --date=short --pretty="%ad" HEAD^..HEAD)
ifeq ("$(RETISTRUCT_GIT_COMMIT_HASH)", "")
	RETISTRUCT_GIT_COMMIT_HASH=1000
endif
RETISTRUCT_PACKAGE=retistruct_$(RETISTRUCT_VERSION).tar.gz
RETISTRUCTDEMOS_PACKAGE=retistructdemos_$(RETISTRUCTDEMOS_VERSION).tar.gz

dist: retistruct check user-guide
	rm -f retistruct_$(RETISTRUCT_VERSION).zip
	cp retistruct.Rcheck/retistruct-manual.pdf doc
	zip -r retistruct_$(RETISTRUCT_VERSION).zip $(RETISTRUCT_PACKAGE) $(RETISTRUCTDEMOS_PACKAGE) matlab  doc/retistruct-user-guide.pdf doc/retistruct-manual.pdf install.R README.md  -x '*/*.svn/*' -x 'matlab/images/*'

roxygen:
	rm -f pkg/retistruct/man/*
	echo "if (!library(roxygen2, logical.return=TRUE)) {install.packages(\"roxygen2\"); library(roxygen2) } ; roxygenize(\"pkg/retistruct\", clean=TRUE)" |	R --no-restore --slave

fix-revision:
	perl -p -i -e "s/^retistruct\.global\.revision.*/retistruct.global.revision <- \"$(RETISTRUCT_GIT_COMMIT_HASH)\"/;" pkg/retistruct/R/revision.R
	perl -p -i -e "s/^Date:.*/Date: $(RETISTRUCT_GIT_AUTHORDATE)/;" pkg/retistruct/DESCRIPTION

retistruct: fix-revision roxygen 
	R CMD build pkg/retistruct

retistructdemos: retistruct
	R CMD build pkg/retistructdemos

install: install-retistruct  install-retistructdemos

install-retistruct: retistruct
	R CMD INSTALL --latex $(RETISTRUCT_PACKAGE) 

install-retistructdemos: retistructdemos
	R CMD INSTALL --latex $(RETISTRUCTDEMOS_PACKAGE) 

doc: retistruct
	rm -f retistruct.pdf
	R CMD Rd2dvi --pdf retistruct

user-guide:
	cd doc &&	pdflatex retistruct-user-guide.tex && cp retistruct-user-guide.pdf ../www

check:
	R CMD check --as-cran $(RETISTRUCT_PACKAGE)

revision:
	@echo $(RETISTRUCT_GIT_COMMIT_HASH)

clean:
	rm -f pkg/retistruct/R/*~
	rm -f trunk/matlab/*~
