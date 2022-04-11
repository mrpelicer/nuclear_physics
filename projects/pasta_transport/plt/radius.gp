set output 'radius.tex'
set terminal cairolatex pdf color colortext standalone size 10 cm, 7 cm

# set multiplot layout 1,2 rowsfirst

set xlabel "$\\rho_B$(fm$^{-3}$)"

# set linetype 1 linecolor rgb "#99004d"
# set linetype 2 linecolor rgb "#cca300"
# set linetype 3 linecolor rgb "#3333cc"
# set linetype 4 linecolor rgb "#cc6600"
# set linetype 5 linecolor rgb "#33cc33"

set xtics 0.02
# set ytics 2
#set xtics 0.02
set xrange[0:.09]
#set yrange[1:9]

set key at 0.035,6.8 spacing 1.3

set ylabel offset 0.0,0.0 "$R_d$\\,(fm)"
set cbrange [1:5]

# define the palette to your liking
set palette defined( 1.0 "#99004d", 2.0 "#cca300", 3.0 "#3333cc", 4.0 "#cc6600", 5.0 "#33cc33")

unset colorbox

set ylabel offset 0.0,0.0 "$R_d$\\, (fm)"
plot 	"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:9:13 every 4 w lp pt 6  ps .8 lw 6 lc palette notitle ,\
	"../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:9:13 every 4 w lp pt 12 ps .9 lw 6 lc palette notitle ,\
       NaN lc rgb "#99004d" lw 8 title "droplets" ,\
    	NaN lc rgb "#cca300" lw 8 title "rods",\
    	NaN lc rgb "#3333cc" lw 8 title "slabs" ,\
    	NaN lc rgb "#cc6600" lw 8 title "tubes",\
    	NaN lc rgb "#33cc33" lw 8 title "bubbles" ,\
       NaN w lp pt 6  ps .8 lw 6 lc rgb "black" title "T=1 MeV",\
       NaN w lp pt 12 ps .9 lw 6 lc rgb "black" title "T=3 MeV"



unset multiplot


