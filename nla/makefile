filename=homg
abstract=homg_abstract
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


# eof
