.SUFFIXES: .tex .odt .pdf

FILE=retistruct_plos_oo

.odt.tex:
	w2l -config plos.xml -inputencoding latin1 -ignore_empty_paragraphs true -wrap_lines_after 70 -no_preamble=true -formulas latex  $<
	cp $@ `basename $< .odt`-raw.tex
	./cleanlatex.pl <  `basename $< .odt`-raw.tex > `basename $< .odt`-clean.tex
	cat plos_header.tex `basename $< .odt`-clean.tex plos_footer.tex > $@

.tex.pdf:
	pdflatex $<
	bibtex `basename $< .tex`
	pdflatex $<
	pdflatex $<
	bibtex `basename $< .tex`

clean: 
	rm -f *.aux

all:
	w2l -config plos.xml -inputencoding latin1 -ignore_empty_paragraphs true -wrap_lines_after 70 -no_preamble=true -formulas latex  $(FILE).odt  ##  
	./cleanlatex.pl <  $(FILE).tex > $(FILE)-clean.tex
	cat plos_header.tex $(FILE)-clean.tex plos_footer.tex > $(FILE)-out.tex	
	pdflatex $(FILE)-out.tex
	bibtex $(FILE)-out
	pdflatex $(FILE)-out.tex
	pdflatex $(FILE)-out.tex
	bibtex $(FILE)-out



