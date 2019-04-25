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

FORCE: build/aufg2-energy.pdf \
	build/aufg2-position.pdf | build

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -pedantic -o build/aufg1

aufg1: build/aufg1 | build
	build/aufg1

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -pedantic -o build/aufg2

aufg2: build/aufg2 | build
	build/aufg2

build/aufg2-energy.pdf build/aufg2-position.pdf: aufgabe2-plot.py build/aufg2-energy.txt build/aufg2-position.txt | build
	python aufgabe2-plot.py

build/aufg2-energy.txt build/aufg2-position.txt: build/aufg2 | build
	build/aufg2

build:
	mkdir -p build/

clean:
	rm -rf build
