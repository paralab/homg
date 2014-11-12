filename1=reply1
filename2=reply2
compile= ${compPdftex}
compPdftex=pdflatex

all:
	${compile}  ${filename1}
	-bibtex ${filename1}
	bibtex ${filename1}
	${compile}  ${filename2}
	-bibtex ${filename2}
	bibtex ${filename2}
	@echo Done ....!

clean:
	rm -f *.aux *.log *.dvi *.blg *.bbl *.out
	@echo Cleaned ....! 

# eof
