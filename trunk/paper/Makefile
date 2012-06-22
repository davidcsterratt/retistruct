all:
	w2l -clean -inputencoding latin1 -ignore_empty_paragraphs true -wrap_lines_after 70 -formulas latex retistruct_plos_oo.odt
	perl -i -p -e 's|\\textit\{(\w+)\}(\w)|\$\1\$\2|g;' retistruct_plos_oo.tex
	pdflatex retistruct_plos_oo.tex

