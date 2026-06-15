all: pdf clean

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

## use texdiff to show difference to previous version
$(FILE)-changes.pdf: tex old/$(FILE)-v2.tex old/ld.cfg
	-latexdiff -c old/ld.cfg old/$(FILE)-v2.tex $(FILE).tex > $(FILE)-changes.tex
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-bibtex $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes

## use texdiff to show difference to previous version
supplement-changes.pdf: supplement old/supplement-v2.tex old/ld.cfg
	-latexdiff -c old/ld.cfg old/supplement-v2.tex supplement.tex > supplement-changes.tex
	-pdflatex --interaction nonstopmode supplement-changes
	-bibtex supplement-changes
	-pdflatex --interaction nonstopmode supplement-changes
	-pdflatex --interaction nonstopmode supplement-changes

## clean
clean:  
	-rm $(FILE).aux  $(FILE).blg  $(FILE).log $(FILE).bbl $(FILE).fff \
	$(FILE).out $(FILE).brf $(FILE).tex $(FILE).synctex.gz \
	supplement.aux  supplement.blg  supplement.log  supplement.tex  supplement.bbl \
	supplement.out supplement.brf supplement.synctex.gz \
	$(FILE)-changes.aux  $(FILE)-changes.blg  $(FILE)-changes.log $(FILE)-changes.bbl $(FILE)-changes.fff \
	$(FILE)-changes.out $(FILE)-changes.brf $(FILE)-changes.tex $(FILE)-changes.synctex.gz \
	supplement-changes.aux  supplement-changes.blg  supplement-changes.log supplement-changes.bbl supplement-changes.fff \
	supplement-changes.out supplement-changes.brf supplement-changes.tex supplement-changes.synctex.gz

cleanall: clean
	-rm -r ./figure/
	-rm -r ./cache/

.PHONY: all clean cleanall pdf code pdf2 tex
