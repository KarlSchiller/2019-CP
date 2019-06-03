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

FORCE: build/aufg1_a1.pdf \
	build/aufg1_a2.pdf \
	build/aufg1_a3.pdf \
	build/aufg1_b.pdf

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -pedantic -o build/aufg1

aufg1: build/aufg1 | build
	build/aufg1

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -pedantic -o build/aufg2

aufg2: build/aufg2 | build
	build/aufg2

build/aufg1_a1.pdf build/aufg1_a2.pdf build/aufg1_a3.pdf build/aufg1_b.pdf: aufgabe1-plot.py build/aufg1_a1.txt build/aufg1_a2.txt build/aufg1_a3.txt build/aufg1_b.txt | build
	python aufgabe1-plot.py

build/aufg1_a1.txt build/aufg1_a2.txt build/aufg1_a3.txt build/aufg1_b.txt: build/aufg1 | build
	build/aufg1

build:
	mkdir -p build/

clean:
	rm -rf build
