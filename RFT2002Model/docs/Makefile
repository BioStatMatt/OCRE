rodriguez4-notes.pdf:	rodriguez4-notes.Rnw
	R --vanilla -e "Sweave('rodriguez4-notes.Rnw')"
	R CMD pdflatex rodriguez4-notes.tex

clean:	
	rm -rf rodriguez4-notes.log
	rm -rf rodriguez4-notes-*.pdf
	rm -rf rodriguez4-notes.out
	rm -rf rodriguez4-notes.tex
	rm -rf rodriguez4-notes.aux
	rm -rf rodriguez4-notes.pdf
	rm -rf rodriguez4-notes-concordance.tex

