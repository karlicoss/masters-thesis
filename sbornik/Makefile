NAME := sbornik_gerasimovda.tex

all:
	latexmk -outdir=build -shell-escape -pdf $(NAME)

preview:
	latexmk -outdir=build -shell-escape -pdf -pvc $(NAME)


clean:
	latexmk -outdir=build -c
