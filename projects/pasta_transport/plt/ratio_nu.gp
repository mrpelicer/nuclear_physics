set output 'ratio_nu.tex'
set terminal cairolatex pdf color colortext standalone size 18 cm, 8 cm

set multiplot layout 1,2 rowsfirst

set yrange[0.8:1.15]
set key top right spacing 1.5

set rmargin 2
set lmargin 11
set title "Rod: $n_B =0.06$ fm$^{-3}$"
set xlabel "$L_2/R_2$"
set arrow 1 from 1,1 to 10000,1 nohead dt 2
set log x
set ylabel offset 2.0,0.0 "\\nu_{p}/\\nu_{a}" rotate by 0
plot 	"../data/transport_Ld_T1.000000_2.txt"   			u 1:5 w l  lw 8 lc rgb "blue" title "1 MeV" ,\
			"../data/transport_Ld_T3.000000_2.txt"   			u 1:5 w l  lw 8 lc rgb "red" 	title "3 MeV"

set title "Slab: $n_B =0.08$ fm$^{-3}$"
unset ylabel
unset rmargin
unset key 
set lmargin 8
set xlabel "$L_1/R_1$"
# set ylabel offset 2.0,0.0 "\\nu_{p}/\\nu_{a}" rotate by 0
plot 	"../data/transport_Ld_T1.000000_1.txt"   			u 1:5 w l lw 8 lc rgb "blue" title "1 MeV" ,\
			"../data/transport_Ld_T3.000000_1.txt"   			u 1:5 w l lw 8 lc rgb "red" 	title "3 MeV"