gnuplot $1.gp
pdflatex --interaction=batchmode $1.tex

rm $1.log $1.tex $1.aux $1-inc.pdf $1.dvi *~

xdg-open $1.pdf
