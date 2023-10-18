all: doc figs

doc: documentation.tex
	pdflatex documentation.tex
	pdflatex documentation.tex

figs: amr-grids.tex amr-grids.tikz
	pdflatex amr-grids.tex

clean:
	rm -f *.aux *.log *.fdb_latexmk *.fls *.out