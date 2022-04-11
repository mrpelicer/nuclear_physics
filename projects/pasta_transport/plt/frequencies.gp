set output 'frequencies.tex'
set terminal cairolatex pdf color colortext standalone size 30 cm, 20 cm

set multiplot layout 2,2
set xlabel "$\\rho_B$(fm$^{-3}$)"
#set linetype 1 linecolor rgb "#99004d"
#set linetype 2 linecolor rgb "#cca300"
#set linetype 3 linecolor rgb "#3333cc"
#set linetype 4 linecolor rgb "#cc6600"
#set linetype 5 linecolor rgb "#33cc33"
#
#set xtics 0.03
#set ytics ?
# set xrange[:0.95]

set key Left left  bottom spacing 1.5
#set label "$Y_p=0.3$" at 0.09,37 rotate by 5
#set label "$Y_p=0.5$" at 0.09,75 rotate by 15
#set label "$Y_p=0.1$" at 0.09,18 


set ylabel offset 6.0,0.0 "$x_d$" rotate by 0	
plot  "../data/coulomb_T0.000000.txt" u 1:2 	w l lw 6 			lc rgb "#black"   	title "x_{3}" ,\
			"../data/coulomb_T0.000000.txt" u 1:3		w l lw 7 			lc rgb "#cca300"   	title "x_{2a}" ,\
			"../data/coulomb_T0.000000.txt" u 1:4		w l lw 8 			lc rgb "#cca300"   	title "x_{2p}" ,\
			"../data/coulomb_T0.000000.txt" u 1:5		w l lw 9 			lc rgb "#cca300"   	title "x_{1a}" ,\
			"../data/coulomb_T0.000000.txt" u 1:6		w l lw 10			lc rgb "#cca300"   	title "x_{1p}"

set ylabel offset 0.0,0.0 "$\\langle \\nu \\rangle$" rotate by 0	
plot  "../data/coulomb_T0.000000.txt" u 1:($8/$9)		w l lw 7 			lc rgb "#cca300"   	title "$\\nu_{a2}/\\nu_{p2}$" ,\
			"../data/coulomb_T0.000000.txt" u 1:($10/$11)		w l lw 8 			lc rgb "#cca300"   	title "$\\nu_{a1}/\\nu_{p1}$"

set ylabel offset 0.0,0.0 "$\\langle \\nu \\rangle$" rotate by 0	
plot  "../data/coulomb_T0.000000.txt" u 1:7 	w l lw 6 			lc rgb "#black"   	title "x_{3}" ,\
			"../data/coulomb_T0.000000.txt" u 1:8		w l lw 7 			lc rgb "#cca300"   	title "x_{2a}" ,\
			"../data/coulomb_T0.000000.txt" u 1:9		w l lw 8 			lc rgb "#cca300"   	title "x_{2p}" ,\
			"../data/coulomb_T0.000000.txt" u 1:10	w l lw 9 			lc rgb "#cca300"   	title "x_{1a}" ,\
			"../data/coulomb_T0.000000.txt" u 1:11	w l lw 10			lc rgb "#cca300"   	title "x_{1p}"
