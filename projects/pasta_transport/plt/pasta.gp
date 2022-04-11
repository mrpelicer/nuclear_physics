set output 'pasta.tex'
set terminal cairolatex pdf color colortext standalone size 20 cm, 14 cm

set multiplot layout 2,2 rowsfirst

set xlabel "$\\rho_B$(fm$^{-3}$)"

set key at 0.09,6 spacing 1.5

#set ytics 20.
set xtics 0.03
set xrange[:.095]
set yrange[-6:]
set cbrange [1:5]

# define the palette to your liking
set palette defined( 1.0 "#99004d", 2.0 "#cca300", 3.0 "#3333cc", 4.0 "#cc6600", 5.0 "#33cc33")

unset colorbox
set rmargin 6

plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:3:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:3:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle ,\
			NaN lc rgb "#99004d" lw 8 title "droplets" ,\
    	NaN lc rgb "#cca300" lw 8 title "rods",\
    	NaN lc rgb "#3333cc" lw 8 title "slabs" ,\
    	NaN lc rgb "#cc6600" lw 8 title "tubes",\
    	NaN lc rgb "#33cc33" lw 8 title "bubbles" ,\
			NaN w lp pt 6  ps .8 lw 6 lc rgb "black" title "T=1 MeV",\
			NaN w lp pt 12 ps .9 lw 6 lc rgb "black" title "T=3 MeV"

unset yrange
set yrange[:0.15]
set ylabel offset 0.0,0.0 "Y_p"
plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:7:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:7:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle

unset yrange
#set yrange[:0.15]
set ylabel offset 0.0,0.0 "A_e, Z_e"
plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:4:13 every 4 w lp pt 6  dt 2 ps .8 lw 6 lc rgb "black" title "$Z_e$" ,\
			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:4:13 every 4 w lp pt 12 dt 2 ps .9 lw 6 lc rgb "black" notitle	,\
		 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:5:13 every 4 w lp pt 6  dt 1 ps .8 lw 6 lc rgb "black" title "$A_e$" ,\
			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:5:13 every 4 w lp pt 12 dt 1 ps .9 lw 6 lc rgb "black" notitle			

unset yrange
#set yrange[:1e-6]
set ylabel offset 0.0,0.0 "R_d"
plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:9:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
			"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:9:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle
