all: build/dokument.pdf

TeXOptions = -lualatex \
				 -interaction=nonstopmode \
				 -halt-on-error \
				 -output-directory=build

build/dokument.pdf: FORCE | build
	latexmk $(TexOptions) dokument.tex

light: dokument.tex
	lualatex -interaction=nonstopmode \
	-halt-on-error \
	-output-directory=build \
	dokument.tex

FORCE: build/leastsquares.pdf | build

build/leastsquares.pdf: aufgabe2-plot.py | build
	python aufgabe2-plot.py

build:
	mkdir -p build/

clean:
	rm -rf build
