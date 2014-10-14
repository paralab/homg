filename=homg
abstract=homg_abstract
filename_diff=${filename}_diff
compilePdfOptions=-shell-escape

compile= ${compPdftex}
compPdftex=pdflatex ${compilePdfOptions}

all:
	${compile}  ${filename}
	-bibtex ${filename}
	bibtex ${filename}
#	${compile} ${filename}
#	${compile} ${filename}
#	${compile} ${filename}
	@echo Done ....!

pdf:
	${compPdftex} ${filename}

abstract:
	${compPdftex} ${abstract}

bib:
	${compile} ${filename}
	bibtex ${filename}

clean:
	rm -f *.aux *.log *.dvi *.blg *.bbl *.out ${filename}-figure*
	@echo Cleaned ....! 

%.pdf: %.tex
	${compile} $*
	${compile} $* 

diff:
	latexdiff ${filename}_submitted.tex ${filename}.tex > ${filename_diff}.tex
	${compile} ${filename_diff}
	-bibtex ${filename_diff}
	bibtex ${filename_diff}

# eof
