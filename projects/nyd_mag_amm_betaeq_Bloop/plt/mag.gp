set output 'mag.tex'
set terminal cairolatex pdf color colortext standalone size 18 cm, 8 cm

#set multiplot layout 2,2 rowsfirst

set xlabel "$B(G)$"
set ylabel offset 0.0,0.0 "$M(G)$" rotate by 90
# set format y "$10^{%T}$"

# set label "L3\\omega\\rho+\\alpha_v=1.0+$\\beta=1.2" at .05,0.6
#set label "B=3 $\\times 10^{18}$ G, SEM momento magn." at .3,0.6
set rmargin 60
# set log y
#set yrange[-1e16:]
set zeroaxis
set xrange[1e17:]
set log x
set rmargin 20
set key left
# set key at 1.42,0.6 Left spacing 2.

parameter="l3wr"

#Plot	1
plot  "../data/magnetization_".parameter."_wtB.txt"   	u 1:3  w l lw 10 dt 2 lc rgb "black"   title "wtB" ,\
			 "../data/magnetization_".parameter."_wtA.txt"   	u 1:3  w l lw 10 dt 1 lc rgb "black"   title "wtA"

			 #"../data/magnetization_".parameter."_wtB.txt"   	u 1:4  w l lw 5 dt 2 lc rgb "#00ff99"   title "wtB",\
			 #"../data/magnetization_".parameter."_wtA.txt"   	u 1:4  w l lw 5 dt 1 lc rgb "#00ff99"   title "wtA"