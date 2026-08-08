.PHONY: all clean

all: main.pdf

main.pdf: main.tex references.bib $(wildcard figs/*.tex)
	lualatex -interaction=nonstopmode main.tex
	bibtex main
	lualatex -interaction=nonstopmode main.tex
	lualatex -interaction=nonstopmode main.tex

clean:
	rm -f *.aux *.log *.out *.bbl *.blg
