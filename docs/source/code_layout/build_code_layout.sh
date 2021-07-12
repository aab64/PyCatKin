#!/bin/bash

file="code_layout.svg"
if [ -f "$file" ] ; then
    rm "$file"
fi

latex code_layout.tex
dvisvgm  --exact --font-format=woff code_layout.dvi
pdflatex code_layout.tex

rm -f *.aux *.gz *.log *.dvi