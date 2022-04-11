set output 'dens_y.tex'
set terminal cairolatex pdf color colortext standalone size 18 cm, 8 cm

#set multiplot layout 2,2 rowsfirst

set xlabel "$\\rho_B$ (fm$^{-3}$)"
set ylabel offset 0.0,0.0 "$\\rm \\rho_i$ (fm$^{-3}$)" rotate by 90
set format y "$10^{%T}$"

set label "L3\\omega\\rho+\\alpha_v=1.0+$\\beta=1.2" at .05,0.6
#set label "B=3 $\\times 10^{18}$ G, SEM momento magn." at .3,0.6
set rmargin 60
set log y
set yrange[1e-4:]
set xrange[0.:1.]

set rmargin 20
set key at 1.25,0.6 Left spacing 2.

parameter="l3wr"
# "../data/dens_beta_".parameter."_noB.txt"  		u 1:2  w l lw 7	 dt 3 lc rgb "#ffcc00"   	title "p" ,\
	 		# "../data/dens_beta_".parameter."_noB.txt"   	u 1:3  w l lw 7	 dt 3 lc rgb "#00ff99"   	title "n" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:4  w l lw 7	 dt 3	lc rgb "#0066cc"   	title "e" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:5  w l lw 7	 dt 3	lc rgb "#660066"   	title "$\\mu$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:6  w l lw 7	 dt 3	lc rgb "#33cccc"   	title "$\\Lambda^0$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:7  w l lw 7	 dt 3	lc rgb "#6600ff"   	title "$\\Sigma^+$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:8  w l lw 7	 dt 3	lc rgb "#bf00ff"   	title "$\\Sigma^0$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:9  w l lw 7	 dt 3	lc rgb "#ff9933"   	title "$\\Sigma^-$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:10 w l lw 7	 dt 3	lc rgb "#00cc00"   	title "$\\Xi^0$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:11 w l lw 7	 dt 3	lc rgb "#ff3300"   	title "$\\Xi^-$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:12 w l lw 7	 dt 3	lc rgb "#993333"   	title "$\\Delta^{++}$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:13 w l lw 7	 dt 3	lc rgb "#ff66cc"   	title "$\\Delta^+$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:14 w l lw 7	 dt 3	lc rgb "#006600"   	title "$\\Delta^0$" ,\
			# "../data/dens_beta_".parameter."_noB.txt"   	u 1:15 w l lw 7	 dt 3	lc rgb "#000000"    title "$\\Delta^-$" ,\

#Plot	1
plot  "../data/dens_beta_".parameter."_wtB.txt"   	u 1:2  w l lw 12 dt 2 lc rgb "#00ff99"   title "p" ,\
			"../data/dens_beta_".parameter."_wtB.txt"  		u 1:3  w l lw 12 dt 2 lc rgb "#ffcc00"   title "n" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:4  w l lw 12 dt 2 lc rgb "#0066cc"   title "e" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:5  w l lw 12 dt 2 lc rgb "#660066"   title "$\\mu$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:6  w l lw 12 dt 2 lc rgb "#33cccc"   title "$\\Lambda^0$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:7  w l lw 12 dt 2 lc rgb "#6600ff"   title "$\\Sigma^+$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:8  w l lw 12 dt 2 lc rgb "#bf00ff"   title "$\\Sigma^0$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:9  w l lw 12 dt 2 lc rgb "#ff9933"   title "$\\Sigma^-$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:10 w l lw 12 dt 2 lc rgb "#00cc00"   title "$\\Xi^0$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:11 w l lw 12 dt 2 lc rgb "#ff3300"   title "$\\Xi^-$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:12 w l lw 12 dt 2 lc rgb "#993333"   title "$\\Delta^{++}$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:13 w l lw 12 dt 2 lc rgb "#ff66cc"   title "$\\Delta^+$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:14 w l lw 12 dt 2 lc rgb "#006600"   title "$\\Delta^0$" ,\
			"../data/dens_beta_".parameter."_wtB.txt"   	u 1:15 w l lw 12 dt 2 lc rgb "#000000"   title "$\\Delta^-$" ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:2  w l lw 10 dt 1 lc rgb "#00ff99"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"  	u 1:3  w l lw 10 dt 1 lc rgb "#ffcc00"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:4  w l lw 10 dt 1 lc rgb "#0066cc"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:5  w l lw 10 dt 1 lc rgb "#660066"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:6  w l lw 10 dt 1 lc rgb "#33cccc"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:7  w l lw 10 dt 1 lc rgb "#6600ff"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:8  w l lw 10 dt 1 lc rgb "#bf00ff"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:9  w l lw 10 dt 1 lc rgb "#ff9933"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:10 w l lw 10 dt 1 lc rgb "#00cc00"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:11 w l lw 10 dt 1 lc rgb "#ff3300"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:12 w l lw 10 dt 1 lc rgb "#993333"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:13 w l lw 10 dt 1 lc rgb "#ff66cc"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:14 w l lw 10 dt 1 lc rgb "#006600"   notitle ,\
			 "../data/dens_beta_".parameter."_wtA.txt"   	u 1:15 w l lw 10 dt 1 lc rgb "#000000"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"  		u 1:2  w l lw 10 dt 2 lc rgb "#ffcc00"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:3  w l lw 10 dt 2 lc rgb "#00ff99"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:4  w l lw 10 dt 2 lc rgb "#0066cc"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:5  w l lw 10 dt 2 lc rgb "#660066"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:6  w l lw 10 dt 2 lc rgb "#33cccc"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:7  w l lw 10 dt 2 lc rgb "#6600ff"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:8  w l lw 10 dt 2 lc rgb "#bf00ff"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:9  w l lw 10 dt 2 lc rgb "#ff9933"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:10 w l lw 10 dt 2 lc rgb "#00cc00"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:11 w l lw 10 dt 2 lc rgb "#ff3300"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:12 w l lw 10 dt 2 lc rgb "#993333"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:13 w l lw 10 dt 2 lc rgb "#ff66cc"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:14 w l lw 10 dt 2 lc rgb "#006600"   notitle ,\
#			"../data/dens_beta_".parameter."_wtB1.txt"   	u 1:15 w l lw 10 dt 2 lc rgb "#000000"   notitldata
