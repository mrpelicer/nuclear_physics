set output 'radius.tex'
set terminal cairolatex pdf color colortext standalone size 15 cm, 7 cm

set multiplot layout 1,2 rowsfirst

set xlabel "$\\rho_B$(fm$^{-3}$)"
set ylabel offset 0.0,0.0 "$R_d$\\,(fm)"

set xtics 0.03
set ytics 2
#set xtics 0.02
set xrange[:0.16]
set yrange[1:9]

set key bottom right Right spacing 1.3
set label "$Y_p=0.5$" at 0.07,70 rotate by 10
set label "$Y_p=0.3$" at 0.06,35 rotate by 5

set title "$T=1$ MeV"

plot "../../data/data_cpa/cpa_iufsu_yp0.500000_T1.000000.txt"   	u 1:10:12 w p ps 0.6 pt 7 lc variable title '$Y_p=0.5$',\
	  "../../data/data_cpa/cpa_iufsu_yp0.300000_T1.000000.txt"   	u 1:10:12 w p ps 0.6 pt 4 lc variable title '$Y_p=0.3$'

set title "$T=3$ MeV"
plot "../../data/data_cpa/cpa_iufsu_yp0.500000_T3.000000.txt"   	u 1:10:12 w p ps 0.6 pt 7 lc variable notitle,\
	  "../../data/data_cpa/cpa_iufsu_yp0.300000_T3.000000.txt"   	u 1:10:12 w p ps 0.6 pt 4 lc variable notitle,\
       NaN lc 1 lw 8 title "spheres" ,\
       NaN lc 2 lw 8 title "rods",\
       NaN lc 3 lw 8 title "slabs" ,\
       NaN lc 4 lw 8 title "tubes",\
       NaN lc 5 lw 8 title "bubbles"


unset multiplot


