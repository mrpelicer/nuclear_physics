set output 'formfactor_L.tex'
set terminal cairolatex pdf color colortext standalone size 32 cm, 13 cm

set multiplot layout 1,3 rowsfirst

set xlabel "q/k_F$"

set key top right spacing 1.5 samplen 10
# set lmargin 8
# set rmargin 0
# set size 0.25,1
# set key rmargin
# set tmargin 3

# # set size ratio -1
id=2

# plot 	"../data/form_factor1.000000_".id.".txt"   			u 1:2 w l dt 3 lw 15 lc rgb "black" title "$ F_a$" ,\
# 			"../data/form_factor1.000000_".id.".txt"   			u 1:3 w l dt 2 lw 15 lc rgb "red" title "$ F_p$"
# set format y ''
# plot 	"../data/form_factor10.000000_".id.".txt"   			u 1:2 w l dt 3 lw 15 lc rgb "black" notitle ,\
# 			"../data/form_factor10.000000_".id.".txt"   			u 1:3 w l dt 2 lw 15 lc rgb "red" notitle

# plot 	"../data/form_factor100.000000_".id.".txt"   			u 1:2 w l dt 3 lw 15 lc rgb "black" notitle ,\
# 			"../data/form_factor100.000000_".id.".txt"   			u 1:3 w l dt 2 lw 15 lc rgb "red" notitle

# plot 	"../data/form_factor1000.000000_".id.".txt"   			u 1:2 w l dt 3 lw 15 lc rgb "black" notitle ,\
# 			"../data/form_factor1000.000000_".id.".txt"   			u 1:3 w l dt 2 lw 15 lc rgb "red" notitle

unset ylabel
set yrange[0:]
set label "\\Large \\langle F_".id."^2\\rangle_a" at 1.3,0.2 
set title " "
plot 	"../data/save/form_factor1.000000_".id.".txt"   			u 1:2 w l dt 1 lw 15 lc rgb "black" 	title "$ L_1= R_1$" ,\
			"../data/save/form_factor10.000000_".id.".txt"   			u 1:2 w l dt 2 lw 15 lc rgb "#800000" 	title "$ L_1= 10 R_1$" ,\
			"../data/save/form_factor100.000000_".id.".txt"   		u 1:2 w l dt 5 lw 15 lc rgb "#006600" title "$ L_1= 100 R_1$" ,\
			"../data/save/form_factor1000.000000_".id.".txt"   		u 1:2 w l dt 3 lw 18 lc rgb "#0033cc" title "$ L_1= 1000 R_1$"

# set format y ''
unset label
set label "\\Large \\langle F_".id."^2\\rangle_p " at 1.3,0.2 

set title "\\Large Rod: $n_B =0.06$ fm$^{-3}$, $T=1$ MeV, $R_2= 4.59$ fm"
# set title "\\Large Slab: $n_B =0.08$ fm$^{-3}$, $T=1$ MeV, $R_1= 3.42$ fm"
plot 	"../data/save/form_factor1.000000_".id.".txt"   			u 1:3 w l dt 1 lw 15 lc rgb "black" 	title "$ L_1= R_1$" ,\
			"../data/save/form_factor10.000000_".id.".txt"   			u 1:3 w l dt 2 lw 15 lc rgb "#800000" 	title "$ L_1= 10 R_1$" ,\
			"../data/save/form_factor100.000000_".id.".txt"   		u 1:3 w l dt 5 lw 15 lc rgb "#006600" title "$ L_1= 100 R_1$" ,\
			"../data/save/form_factor1000.000000_".id.".txt"   		u 1:3 w l dt 3 lw 18 lc rgb "#0033cc" title "$ L_1= 1000 R_1$"

unset label
set label " \\Large \\langle F_".id."^2\\rangle_a +2\\langle F_".id."^2\\rangle_p" at .9,0.57
set title " "
plot 	"../data/save/form_factor1.000000_".id.".txt"   			u 1:4 w l dt 1 lw 15 lc rgb "black" 	title "$ L_1= R_1$" ,\
			"../data/save/form_factor10.000000_".id.".txt"   			u 1:4 w l dt 2 lw 15 lc rgb "#800000" 	title "$ L_1= 10 R_1$" ,\
			"../data/save/form_factor100.000000_".id.".txt"   		u 1:4 w l dt 5 lw 15 lc rgb "#006600" title "$ L_1= 100 R_1$" ,\
			"../data/save/form_factor1000.000000_".id.".txt"   		u 1:4 w l dt 3 lw 18 lc rgb "#0033cc" title "$ L_1= 1000 R_1$"


