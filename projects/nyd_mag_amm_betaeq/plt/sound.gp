set output 'sound.tex'
set terminal cairolatex pdf color colortext standalone size 18 cm, 20 cm

set multiplot layout 4,2 rowsfirst

set xlabel "$\\rho_B$ (fm$^{-3}$)"
#set format y "$10^{%T}$"

#set label "L3\\omega\\rho+\\alpha_v=0.5+$\\beta=1.0" at .05,0.6
#set label "B=3 $\\times 10^{18}$ G, SEM momento magn." at .3,0.6
#set rmargin 60
#set log y
#set yrange[0:]
set xrange[0.:1.]

#set rmargin 20
#set key at 1.42,0.6 Left spacing 2.

parameter="gm1"

#press_thermodynamic
set key Left left
set ylabel offset 0.0,0.0 "$P= -\\varepsilon + \\mu_B \\rho_B" rotate by 90
plot  "../data/hmg_beta_".parameter."_noB.txt"  	u 1:2  w l lw 14	 dt 3 lc rgb "black"  title "$B=0$" ,\
			"../data/hmg_beta_".parameter."_wtB.txt"   	u 1:2  w l lw 11 dt 2 lc rgb "blue"   title "$B\\neq 0 \\; \\kappa_b=0$" ,\
 		  "../data/hmg_beta_".parameter."_wtA.txt"   	u 1:2  w l lw 8 dt 1 lc rgb "red"    title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

#unset key
#set key Left left

set ylabel offset 0.0,0.0 "$$c_S^2$" rotate by 90
plot  "../data/bulk_".parameter."_noB.txt"  	u 1:7  w l lw 14	 dt 3 lc rgb "black"  title "$B=0$" ,\
			"../data/bulk_".parameter."_wtB.txt"   	u 1:7  w l lw 11 dt 2 lc rgb "blue"   	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
 		  "../data/bulk_".parameter."_wtA.txt"   	u 1:7  w l lw 8 dt 1 lc rgb "red"     	title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

set ylabel offset 0.0,0.0 "Mev\\; fm$^{-3}$" rotate by 90
plot  "../data/hmg_beta_".parameter."_wtB.txt"   	u 1:12  w l lw 14 dt 2 lc rgb "blue"    title "$P_{ \\parallel}\\; \\kappa_b=0$" ,\
 		  "../data/hmg_beta_".parameter."_wtB.txt"   	u 1:13  w l lw 11 dt 2 lc rgb "red"    title "$P_{ \\perp} \\; \\kappa_b=0$" ,\
			"../data/hmg_beta_".parameter."_wtA.txt"   	u 1:12  w l lw 14 dt 1 lc rgb "blue"    title "$P_{ \\parallel}\\; \\kappa_b\\neq 0$" ,\
			"../data/hmg_beta_".parameter."_wtA.txt"   	u 1:13  w l lw 11 dt 1 lc rgb "red"    title "$P_{\\perp}\\; \\kappa_b\\neq 0$"

set key Left right bottom
set yrange[-12:12]
set ylabel offset 0.0,0.0 "P_\\perp/P_\\parallel$" rotate by 90
plot  "../data/hmg_beta_".parameter."_wtB.txt"   	u 1:($13/$12)  w l lw 10 dt 2 lc rgb "blue"	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
			"../data/hmg_beta_".parameter."_wtA.txt"   	u 1:($13/$12)  w l lw 8 dt 1 lc rgb "red"    title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

unset yrange
set ylabel offset 0.0,0.0 "c_{S \\parallel}^2" rotate by 90
plot  "../data/bulk_".parameter."_wtB.txt"   	u 1:8  w l lw 11 dt 2 lc rgb "blue"   	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
 		  "../data/bulk_".parameter."_wtA.txt"   	u 1:8  w l lw 8 dt 1 lc rgb "red"     	title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

set ylabel offset 0.0,0.0 "c_{S \\perp}^2" rotate by 90
plot  "../data/bulk_".parameter."_wtB.txt"   	u 1:9  w l lw 11 dt 2 lc rgb "blue"   	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
 		  "../data/bulk_".parameter."_wtA.txt"   	u 1:9  w l lw 8 dt 1 lc rgb "red"     	title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

set ylabel offset 0.0,0.0 "(c_{S\\parallel}^2 + 2.c_{S \\perp}^2)/3" rotate by 90
plot  "../data/bulk_".parameter."_wtB.txt"   	u 1:10  w l lw 11 dt 2 lc rgb "blue"   	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
 		  "../data/bulk_".parameter."_wtA.txt"   	u 1:10  w l lw 8 dt 1 lc rgb "red"     	title "$B\\neq 0 \\; \\kappa_b\\neq 0$"

set ylabel offset 0.0,0.0 "M (G)$" rotate by 90
plot  "../data/hmg_beta_".parameter."_wtB.txt"   	u 1:11  w p pt 1 lc rgb "blue"	title "$B\\neq 0 \\; \\kappa_b=0$" ,\
			
#			"../data/hmg_beta_".parameter."_wtA.txt"   	u 1:11  w p pt 1 lc rgb "red"    title "$B\\neq 0 \\; \\kappa_b\\neq 0$"
