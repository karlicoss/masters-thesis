NAME      := demo.tex
BUILD_DIR := build


all: build

.PHONY: build # latexmk does the job for us
build: 
	latexmk -outdir=$(BUILD_DIR) -shell-escape -pdf $(NAME)


.PHONY: preview # latexmk does the job for us
preview:
	latexmk -outdir=$(BUILD_DIR) -shell-escape -pdf -pvc $(NAME)


.PHONY: clean
clean:
	latexmk -outdir=$(BUILD_DIR) -c


.PHONY: wipe
wipe:
	latexmk -outdir=$(BUILD_DIR) -C


.PHONY: export
export: build
	zip demo.zip $(NAME) build/demo.pdf