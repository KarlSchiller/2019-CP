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

FORCE: build/a_harm.pdf \
	build/a_unharm.pdf

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -pedantic -o build/aufg1

aufg1: build/aufg1 | build
	build/aufg1

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -pedantic -o build/aufg2

aufg2: build/aufg2 | build
	build/aufg2

build/a_harm.pdf build/a_unharm.pdf: aufgabe1-plot.py build/a_harm.txt build/a_unharm.txt | build
	python aufgabe1-plot.py

build/a_harm.txt build/a_unharm.txt: build/aufg1 | build
	build/aufg1

build:
	mkdir -p build/

clean:
	rm -rf build
