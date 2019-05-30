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
	build/a_unharm.pdf \
	build/aufg1_c_energie.jpg \
	build/aufg2_a_ellipse.jpg \
	build/aufg2_a_schmetterling.jpg \
	build/aufg2_b_energie.jpg \
	build/aufg2_b_drehimpuls.jpg \
	build/aufg2_c_09.jpg \
	build/aufg2_c_11.jpg

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

build/aufg1_c_energie.txt: build/aufg1 | build
	build/aufg1

build/aufg1_c_energie.jpg: build/aufg1_c_energie.txt aufgabe1-plot.py
	python aufgabe1-plot.py

build/aufg2_a_ellipse.txt build/aufg2_a_schmetterling.txt: build/aufg2 | build
	build/aufg2

build/aufg2_b_energie.txt build/aufg2_b_drehimpuls.txt: build/aufg2 | build
	build/aufg2

build/aufg2_c_11.txt build/aufg2_c_09.txt: build/aufg2 | build
	build/aufg2

build/aufg2_a_ellipse.jpg build/aufg2_a_schmetterling.jpg :build/aufg2_a_ellipse.txt build/aufg2_a_schmetterling.txt aufgabe2-plot.py | build
	python aufgabe2-plot.py

build/aufg2_b_energie.jpg build/aufg2_b_drehimpuls.jpg: build/aufg2_b_energie.txt build/aufg2_b_drehimpuls.txt aufgabe2-plot.py | build
	python aufgabe2-plot.py

build/aufg2_c_11.jpg build/aufg2_c_09.jpg: build/aufg2_c_11.txt build/aufg2_c_09.txt aufgabe2-plot.py | build
	python aufgabe2-plot.py

build:
	mkdir -p build/

clean:
	rm -rf build
