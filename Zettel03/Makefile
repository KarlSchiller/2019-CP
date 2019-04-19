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

FORCE:

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -o -pedantic build/aufg1

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -o -pedantic build/aufg2

build:
	mkdir -p build/

clean:
	rm -rf build
