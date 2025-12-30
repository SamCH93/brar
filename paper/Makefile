all: pdf code clean

# names of .Rnw .bib. and .bst files
FILE=BFBRAR
BIB=bibliography

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate .pdf file (without first knitting tex)
pdf2:
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

supplement: supplement.Rnw appendix.Rnw tex
	Rscript -e "knitr::knit2pdf('supplement.Rnw')" --vanilla

## extract R code
code: $(FILE).Rnw
	Rscript -e "knitr::purl('$(FILE).Rnw')" --vanilla

## clean
clean:  
	-rm $(FILE).aux  $(FILE).blg  $(FILE).log $(FILE).bbl \
	$(FILE).out $(FILE).brf $(FILE).tex $(FILE).synctex.gz \
	supplement.aux  supplement.blg  supplement.log  supplement.tex  supplement.bbl \
	supplement.out supplement.brf supplement.synctex.gz

cleanall: clean
	-rm -r ./figure/
	-rm -r ./cache/
