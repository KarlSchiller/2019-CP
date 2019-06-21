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

FORCE: build/endkonfig.jpg \
	build/equilibration_energy.jpg \
	build/equilibration_temp.jpg \
	build/equilibration_vSx.jpg \
	build/equilibration_vSy.jpg \
	build/init.jpg \
	build/messung_T1_temp.jpg \
	build/messung_T1e-2_temp.jpg \
	build/messung_T1e2_temp.jpg \
	build/thermostat_equi_energy.jpg

build/aufg1: aufgabe1.cpp | build
	g++ aufgabe1.cpp -pedantic -o build/aufg1

aufg1: build/aufg1 | build
	build/aufg1

build/equilibration_energy.jpg build/equilibration_temp.jpg build/equilibration_vSx.jpg build/equilibration_vSy.jpg: build/equilibration.txt aufgabe1-plot.py | build
	python aufgabe1-plot.py

build/endkonfig.jpg build/init.jpg: build/init.txt build/endkonfig.txt aufgabe1-plot.py | build
	python aufgabe1-plot.py

build/messung_T1_temp.jpg build/messung_T1e2_temp.jpg build/messung_T1e-2_temp.jpg: build/messung_T1.txt build/messung_T1e2.txt build/messung_T1e-2.txt aufgabe1-plot.py | build
	python aufgabe1-plot.py

build/endkonfig.txt build/init.txt build/equilibration.txt build/messung_T1.txt build/messung_T1e2.txt build/messung_T1e-2.txt build/thermostat_equi.txt build/thermostat_mess.txt: build/aufg1 | build
	build/aufg1

build:
	mkdir -p build/

clean:
	rm -rf build
