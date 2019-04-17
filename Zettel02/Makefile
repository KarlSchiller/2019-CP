all: build/dokument.pdf

TeXOptions = -lualatex \
				 -interaction=nonstopmode \
				 -halt-on-error \
				 -output-directory=build

build/dokument.pdf: FORCE | build
	python aufgabe1-plot.py
	latexmk $(TexOptions) dokument.tex

light: dokument.tex
	lualatex -interaction=nonstopmode \
	-halt-on-error \
	-output-directory=build \
	dokument.tex

FORCE: build/aufg2-eigenvalues.pdf \
	build/aufg2-testpic.pdf \
	build/aufg2-k200.pdf \
	build/aufg2-k300.pdf

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

build/aufg2: aufgabe2.cpp | build
	g++ aufgabe2.cpp -o build/aufg2

aufg2: build/aufg2-k200.pdf | build


build:
	mkdir -p build/

clean:
	rm -rf build
