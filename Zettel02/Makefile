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

FORCE: build/aufg2-eigenvalues.pdf \
	build/aufg2-testpic.pdf \
	build/aufg2-k200.pdf \
	build/aufg2-k300.pdf \
	build/10.pdf \
	build/20.pdf \
	build/50.pdf

build/aufg2-eigenvalues.pdf: build/aufg2-eigenvalues.txt | build
	python aufgabe2-plot.py

build/aufg2-k200.pdf: build/aufg2-k200.txt | build
	python aufgabe2-plot.py

build/aufg2-k300.pdf: build/aufg2-k300.txt | build
	python aufgabe2-plot.py

build/aufg2-eigenvalues.txt: build/aufg2 | build
	build/aufg2

build/aufg2-k200.txt: build/aufg2 | build
	build/aufg2

build/aufg2-k300.txt: build/aufg2 | build
	build/aufg2

build/aufg2: aufgabe2.cpp Dateien/service.cpp Dateien/Testdata Dateien/Training | build
	g++ aufgabe2.cpp -o build/aufg2

# just a shortcut for testing aufgabe2.cpp
aufg2: build/aufg2 | build
	build/aufg2

build/10.pdf build/20.pdf build/50.pdf: build/bild_10.txt build/bild_20.txt build/bild_50.txt aufgabe1-plot.py | build
	python aufgabe1-plot.py

build/bild_10.txt build/bild_20.txt build/bild_50.txt: build/aufg1 | build
	build/aufg1

build/aufg1: aufgabe1.cpp Dateien/Bild | build
	g++ aufgabe1.cpp -o build/aufg1

build:
	mkdir -p build/

clean:
	rm -rf build
