#!/bin/sh
#R -e 'library(knitr);library(markdown);knit("mexico-analysis.Rmd"); markdownToHTML("mexico-analysis.md","mexico-analysis.html")'
#
#pandoc -s -i mexico-analysis.md -o mexico-analysis.tex -t latex -S



R -e 'library(knitr);knit("mexico-analysis.Rmd")'

pandoc -s -i mexico-analysis.md -o mexico-analysis.tex -t latex -S
