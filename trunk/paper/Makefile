FILE=retistruct_plos_oo

all:
	w2l -clean -inputencoding latin1 -ignore_empty_paragraphs true -wrap_lines_after 70 -formulas latex -no_preamble=true $(FILE).odt  ##  
	./cleanlatex.pl <  $(FILE).tex > $(FILE)-clean.tex
	cat plos_header.tex $(FILE)-clean.tex plos_footer.tex > $(FILE)-out.tex	
	pdflatex $(FILE)-out.tex
	bibtex $(FILE)-out
	pdflatex $(FILE)-out.tex
	pdflatex $(FILE)-out.tex

