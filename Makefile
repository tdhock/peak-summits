paper.pdf: paper.tex refs.bib figure-Mono27ac-label-error.png
	rm -f *.aux *.bbl
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper
figure-Mono27ac-label-error.png: figure-Mono27ac-label-error.R
	R --no-save < $<
