paper.pdf: paper.tex refs.bib figure-Mono27ac-label-error.png figure-Mono27ac-summits-vs-peaks.png figure-AR1-multimodal.png figure-neuro-training.png
	rm -f *.aux *.bbl
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper
figure-Mono27ac-label-error.png: figure-Mono27ac-label-error.R
	R --no-save < $<
figure-Mono27ac-summits-vs-peaks.png: figure-Mono27ac-summits-vs-peaks.R
	R --no-save < $<
figure-AR1-multimodal.png: figure_7/fig_code/fig7.R
	cd figure_7/fig_code && R --no-save < fig7.R
figure-neuro-training.png: figure_7/fig_code/figure-training.R
	cd figure_7/fig_code && R --no-save < figure-training.R
