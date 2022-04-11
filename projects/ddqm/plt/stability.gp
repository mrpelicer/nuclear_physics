set output 'stability.tex'
set terminal cairolatex pdf color colortext standalone size 20 cm, 14 cm

set xlabel offset 0.0,0.0 "$C$"
set ylabel offset 0.0,0.0 "$\\sqrt D$   (MeV)"

set xtics 0.2
set ytics 5
set linetype 1 linecolor rgb "#660066"
set linetype 2 linecolor rgb "#e65c00"
set linetype 3 linecolor rgb "#00ccff"
set linetype 4 linecolor rgb "#339933"
set linetype 5 linecolor rgb "#ff0000"

#set xrange[:0.13]
#set yrange[2:10]

set key Left box opaque at 0.2,195 spacing 1.5 font ",13"
set border back

#set rmargin 20

plot "../data/stability_ms100.000000_T50.000000.txt"   	u 1:2:3 w p pt 7 ps 1 lc variable notitle ,\
  	   NaN lt 1 lw 8 title "Stable 2F" ,\
       NaN lt 2 lw 8 title "Stable SQM",\
       NaN lt 3 lw 8 title "MetaStable SQM" ,\
       NaN lt 4 lw 8 title "Unstable SQM" ,\
       NaN lt 5 lw 8 title "$m_{up}<0 $ before 1.5 fm $^{-3}$"

