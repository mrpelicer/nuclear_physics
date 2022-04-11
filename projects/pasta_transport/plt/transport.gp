set output 'transport.tex'
set terminal cairolatex pdf color colortext standalone size 30 cm, 20 cm

# set multiplot layout 2,2 rowsfirst

set xlabel "$\\rho_B$(fm$^{-3}$)"
#set linetype 1 linecolor rgb "#99004d"
#set linetype 2 linecolor rgb "#cca300"
#set linetype 3 linecolor rgb "#3333cc"
#set linetype 4 linecolor rgb "#cc6600"
#set linetype 5 linecolor rgb "#33cc33"
#
#set xtics 0.03
#set ytics ?
set xrange[:0.07]

set key Left left  bottom spacing 1.5
#set label "$Y_p=0.3$" at 0.09,37 rotate by 5
#set label "$Y_p=0.5$" at 0.09,75 rotate by 15
#set label "$Y_p=0.1$" at 0.09,18 


set title "$T=0$ MeV"

set ylabel offset 6.0,0.0 "\\Lambda_{ei}" rotate by 0	
plot  "../data/coulomb.txt" u 1:10 	w l lw 20				lc rgb "#black"   	title "gsl" ,\
			"../data/coulomb.txt" u 1:11	w l lw 10 			lc rgb "#cca300"   	title "cuba"

# unset key
# set ylabel offset 6.0,0.0 "$\\frac{1}{\\tau_{\\sigma , \\kappa}} (s^{-1})$" rotate by 0	
# set log y
# plot  "../data/coulomb.txt" u 1:9 		w l lw 10				lc rgb "#99004d"   	title "droplet" ,\
# 			"../data/coulomb.txt" u 1:10		w l lw 8 dt 1 	lc rgb "#cca300"   	title "rod" ,\
# 			"../data/coulomb.txt" u 1:11  	w l lw 8 dt 1 	lc rgb "#3333cc"   	title "slab" ,\
# 			"../data/coulomb.txt" u 1:12		w l lw 8 dt 2 	lc rgb "#cca300"   	title "rod_p" ,\
# 			"../data/coulomb.txt" u 1:13		w l lw 8 dt 3 	lc rgb "#cca300"   	title "rod_a" ,\
# 			"../data/coulomb.txt" u 1:14  	w l lw 8 dt 2 	lc rgb "#3333cc"   	title "slab_p" ,\
# 			"../data/coulomb.txt" u 1:15 		w l lw 8 dt 3 	lc rgb "#3333cc"   	title "slab_a"

# set ylabel offset 6.0,0.0 "$\\sigma \\; (s^{-1})$" rotate by 0	
# unset yrange

# plot  "../data/coulomb.txt" u 1:16 		w l lw 10				lc rgb "#99004d"   	notitle ,\
# 			"../data/coulomb.txt" u 1:17		w l lw 8 dt 1 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:18  	w l lw 8 dt 1 	lc rgb "#3333cc"   	notitle ,\
# 			"../data/coulomb.txt" u 1:19		w l lw 8 dt 2 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:20		w l lw 8 dt 3 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:21  	w l lw 8 dt 2 	lc rgb "#3333cc"   	notitle ,\
# 			"../data/coulomb.txt" u 1:22 		w l lw 8 dt 3 	lc rgb "#3333cc"   	notitle

# set ylabel offset 14.0,0.0 "$\\kappa \\; (\\text{erg K^{-1} cm^{-1} s^{-1})}$" rotate by 0	

# plot  "../data/coulomb.txt" u 1:16 		w l lw 10				lc rgb "#99004d"   	notitle ,\
# 			"../data/coulomb.txt" u 1:17		w l lw 8 dt 1 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:18  	w l lw 8 dt 1 	lc rgb "#3333cc"   	notitle ,\
# 			"../data/coulomb.txt" u 1:19		w l lw 8 dt 2 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:20		w l lw 8 dt 3 	lc rgb "#cca300"   	notitle ,\
# 			"../data/coulomb.txt" u 1:21  	w l lw 8 dt 2 	lc rgb "#3333cc"   	notitle ,\
# 			"../data/coulomb.txt" u 1:22 		w l lw 8 dt 3 	lc rgb "#3333cc"   	notitle