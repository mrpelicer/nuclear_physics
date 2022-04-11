set output 'pasta_size.tex'
set terminal cairolatex pdf color colortext standalone size 20 cm, 14 cm

set multiplot layout 2,2 rowsfirst

set xlabel "$L_d/R_d$"

set key top left Left spacing 1.5
set rmargin 6
# set log y
# set log x
set ylabel offset 0.0,0.0 "s^{-1}" rotate by 0
set log x
set log y
plot 	"../data/transport_Ld_T1.000000.txt"   			u 1:2 w l lw 4 lc rgb "black" title "$\\nu_{3d}$",\
			"../data/transport_Ld_T1.000000.txt"   			u 1:3 w lp pt 6  ps .8 lw 6 lc rgb "blue" title "$\\nu_a$" ,\
			"../data/transport_Ld_T1.000000.txt"   			u 1:4 w lp pt 6  ps .8 lw 6 lc rgb "red" 	title "$\\nu_p",\
			"../data/transport_Ld_T1.000000.txt"   			u 1:6 w lp pt 6  ps .8 lw 6 lc rgb "black" 	title "$\\langle \\nu \\rangle"

set ylabel offset 0.0,0.0 "$Z" rotate by 0
plot 	"../data/transport_Ld_T1.000000.txt"   			u 1:13 w lp pt 6  ps .8 lw 6 lc rgb "blue" title "1 MeV"

unset log y
unset ylabel
plot 	"../data/transport_Ld_T1.000000.txt"   			u 1:9   w lp pt 6  ps .8 lw 6 lc rgb "blue" title "$\\Lambda_{C, a}$",\
			"../data/transport_Ld_T1.000000.txt"   			u 1:10  w lp pt 12  ps .8 lw 6 lc rgb "red" title "$\\Lambda_{C, p}$" ,\

set ylabel offset 0.0,0.0 "\\nu_p/\\nu_a" rotate by 0
plot 	"../data/transport_Ld_T1.000000.txt"   			u 1:5 w lp pt 6  ps .8 lw 6 lc rgb "blue" title "1 MeV" ,\
			"../data/transport_Ld_T3.000000.txt"   			u 1:5 w lp pt 12  ps .8 lw 6 lc rgb "red" 	title "3 MeV"

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
