NAME := sbornik_gerasimovda.tex

build:
	latexmk -outdir=build -shell-escape -pdf $(NAME)

preview:
	latexmk -outdir=build -shell-escape -pdf -pvc $(NAME)


clean:
	latexmk -outdir=build -c


export: build
	zip paper.zip graph.eps nsart_eng.cls sbornik_gerasimovda.tex build/sbornik_gerasimovda.pdf
