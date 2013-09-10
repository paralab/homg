filename=homg
abstract=homg_abstract
compilePdfOptions=

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
	rm -f *.aux *.log *.dvi *.blg *.bbl *.out
	@echo Cleaned ....! 

%.pdf: %.tex
	${compile} $*
	${compile} $* 


# eof
