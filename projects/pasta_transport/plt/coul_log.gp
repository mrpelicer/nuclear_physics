set output 'coul_log.tex'
set terminal cairolatex pdf color colortext standalone size 18 cm, 8 cm

set multiplot layout 1,2 rowsfirst

set xlabel "$L_2/R_2$"
set yrange[0.2:0.6]
set ytics 0.1
set key top right spacing 1.5
# set rmargin 6
set title "Rods, $n_B =0.06$ fm$^{-3}$"
set ylabel "\\Lambda_c" offset "0.0,0.0" rotate by 0
set log x
plot 	"../data/transport_Ld_T1.000000_2.txt"   			u 1:9   w l dt 1 lw 6 lc rgb "red"  notitle,\
			"../data/transport_Ld_T1.000000_2.txt"   			u 1:10  w l dt 2 lw 6 lc rgb "red"  notitle,\
			"../data/transport_Ld_T3.000000_2.txt"   			u 1:9   w l dt 1 lw 6 lc rgb "blue" notitle,\
			"../data/transport_Ld_T3.000000_2.txt"   			u 1:10  w l dt 2 lw 6 lc rgb "blue" notitle,\
			NaN  w l dt 1 lw 6 lc rgb "red"  title "T=1 MeV",\
			NaN  w l dt 1 lw 6 lc rgb "blue" title "T=3 MeV" ,\
			NaN  w l dt 1 lw 6 lc rgb "white"  title " ",\
			NaN  w l dt 1 lw 6 lc rgb "black" title  "$\\Lambda_{C \\, a}$",\
			NaN  w l dt 2 lw 6 lc rgb "black" title  "$\\Lambda_{C \\, p}$"

set yrange[0.1:0.7]
set title "Slabs, $n_B =0.08$ fm$^{-3}$"
unset ylabel
set xlabel "$L_1/R_1$"

plot 	"../data/transport_Ld_T1.000000_1.txt"   			u 1:9   w l dt 1 lw 6 lc rgb "red"  notitle,\
			"../data/transport_Ld_T1.000000_1.txt"   			u 1:10  w l dt 2 lw 6 lc rgb "red"  notitle,\
			"../data/transport_Ld_T3.000000_1.txt"   			u 1:9   w l dt 1 lw 6 lc rgb "blue" notitle,\
			"../data/transport_Ld_T3.000000_1.txt"   			u 1:10  w l dt 2 lw 6 lc rgb "blue" notitle,\
			NaN  w l dt 1 lw 6 lc rgb "red"  title "T=1 MeV",\
			NaN  w l dt 1 lw 6 lc rgb "blue" title "T=3 MeV" ,\
			NaN  w l dt 1 lw 6 lc rgb "white"  title " ",\
			NaN  w l dt 1 lw 6 lc rgb "black" title  "$\\Lambda_{C \\, a}$",\
			NaN  w l dt 2 lw 6 lc rgb "black" title  "$\\Lambda_{C \\, p}$"

# set arrow 1 from 1,1 to 10000,1 nohead dt 2

# set ylabel offset 0.0,0.0 "\\nu_{p}/\\nu_{a}" rotate by 0
# plot 	"../data/transport_Ld_T1.000000.txt"   			u 1:5 w lp pt 6  ps .8 lw 6 lc rgb "blue" title "1 MeV" ,\
# 			"../data/transport_Ld_T3.000000.txt"   			u 1:5 w lp pt 12  ps .8 lw 6 lc rgb "red" 	title "3 MeV"

# unset yrange
# set yrange[:0.15]
# set ylabel offset 0.0,0.0 "Y_p"
# plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:7:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
# 			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:7:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle

# unset yrange
# #set yrange[:0.15]
# set ylabel offset 0.0,0.0 "A_e, Z_e"
# plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:4:13 every 4 w lp pt 6  dt 2 ps .8 lw 6 lc rgb "black" title "$Z_e$" ,\
# 			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:4:13 every 4 w lp pt 12 dt 2 ps .9 lw 6 lc rgb "black" notitle	,\
# 		 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:5:13 every 4 w lp pt 6  dt 1 ps .8 lw 6 lc rgb "black" title "$A_e$" ,\
# 			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:5:13 every 4 w lp pt 12 dt 1 ps .9 lw 6 lc rgb "black" notitle			

# unset yrange
# #set yrange[:1e-6]
# set ylabel offset 0.0,0.0 "R_d"
# plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:9:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
# 			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:9:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle
