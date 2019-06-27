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

FORCE: build/aufg2.pdf \
	build/init.pdf \
	build/orte.pdf \
	build/d9S1e2.pdf \
	build/d9S1e4.pdf \
	build/d99S1e3.pdf \
	build/d999S1e2.pdf

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -pedantic -o build/aufg1

aufg1: build/aufg1 | build
	build/aufg1

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -pedantic -o build/aufg2

aufg2: build/aufg2 | build
	build/aufg2

build/aufg3: aufgabe3.cpp | build
	g++ aufgabe3.cpp -pedantic -o build/aufg3

aufg3: build/aufg3 | build
	build/aufg3

build/d9S1e2.pdf build/d9S1e4.pdf build/d99S1e3.pdf build/d999S1e2.pdf: aufgabe1-plot.py build/d9S1e2.txt build/d9S1e4.txt build/d99S1e3.txt build/d999S1e2.txt | build
	python aufgabe1-plot.py

build/d9S1e2.txt build/d9S1e4.txt build/d99S1e3.txt build/d999S1e2.txt: build/aufg1 | build
	build/aufg1

build/init.pdf build/orte.pdf: build/init.txt build/orte.txt aufgabe1-plot.py | build
	python aufgabe1-plot.py

build/init.txt build/orte.txt: build/aufg1 | build
	build/aufg1

build/aufg2.pdf build/aufg2_100.pdf: build/aufg2.txt build/aufg2_100.txt aufgabe2-plot.py | build
	python aufgabe2-plot.py

build/aufg2.txt build/aufg2_100.txt: build/aufg2 | build
	build/aufg2

build:
	mkdir -p build/

clean:
	rm -rf build
