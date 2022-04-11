set output 'geometries.tex'
set terminal cairolatex pdf color colortext standalone size 20 cm, 13 cm

set multiplot layout 2,2 rowsfirst

set xlabel "$\\rho_B$(fm$^{-3}$)"
set ylabel offset 0.0,0.0 "${\\mathcal F}/\\rho_B - M_N$\\,(MeV)"

set linetype 1 linecolor rgb "#99004d"
set linetype 2 linecolor rgb "#cca300"
set linetype 3 linecolor rgb "#3333cc"
set linetype 4 linecolor rgb "#cc6600"
set linetype 5 linecolor rgb "#33cc33"

set xtics 0.03
set xrange[:0.15]
set yrange[-10:100]

set key at 0.02,90 spacing 1.5
set label "$Y_p=0.5$" at 0.09,75 rotate by 15
set label "$Y_p=0.3$" at 0.09,37 rotate by 5
set label "$Y_p=0.1$" at 0.09,18 

set ytics 20.


set title "$T=0$ MeV"
set rmargin 6
plot 	"../data/cpa_iufsu_yp0.500000_T0.000000.txt"   	u 1:3:12 w l lw 7 lc variable notitle,\
	  	"../data/cpa_iufsu_yp0.300000_T0.000000.txt"   	u 1:3:12 w l lw 7 lc variable notitle ,\
	  	"../data/cpa_iufsu_yp0.100000_T0.000000.txt"   	u 1:3:12 w l lw 7 lc variable notitle ,\
			"../data/cpa_iufsu_betaEq_T0.000000.txt"   			u 1:3:12 w l lw 10 lc variable notitle ,\
			NaN lc 1 lw 8 title "droplets" ,\
    	NaN lc 2 lw 8 title "rods",\
    	NaN lc 3 lw 8 title "slabs" ,\
    	NaN lc 4 lw 8 title "tubes",\
    	NaN lc 5 lw 8 title "bubbles"

set title "$T=1$ MeV"
set rmargin 6
plot  "../data/cpa_iufsu_yp0.500000_T1.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle,\
	    "../data/cpa_iufsu_yp0.300000_T1.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
	    "../data/cpa_iufsu_yp0.100000_T1.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
			"../data/cpa_iufsu_betaEq_T1.000000.txt"   			u 1:3:12 w l lw 10 lc variable notitle ,\

#set format y ''
#unset ylabel
#set lmargin 9
#set rmargin 2

set title "$T=3$ MeV"
plot "../data/cpa_iufsu_yp0.500000_T3.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
	 	 "../data/cpa_iufsu_yp0.300000_T3.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
	   "../data/cpa_iufsu_yp0.100000_T3.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
		 "../data/cpa_iufsu_betaEq_T3.000000.txt"   			u 1:3:12 w l lw 10 lc variable notitle ,\


set title "$T=5$ MeV"
plot "../data/cpa_iufsu_yp0.500000_T5.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
	   "../data/cpa_iufsu_yp0.300000_T5.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
	   "../data/cpa_iufsu_yp0.100000_T5.000000.txt"   	u 1:3:12 w l lw 10 lc variable notitle ,\
		 "../data/cpa_iufsu_betaEq_T5.000000.txt"   			u 1:3:12 w l lw 10 lc variable notitle ,\

