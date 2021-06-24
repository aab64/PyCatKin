#!/bin/bash

file="code_layout.svg"
if [ -f "$file" ] ; then
    rm "$file"
fi

# lualatex --output-format=dvi code_layout.tex
latex code_layout.tex
dvisvgm code_layout.dvi

rm -f *.aux *.gz *.log *.dvi