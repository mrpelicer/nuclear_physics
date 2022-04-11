#set output 'form_factors.tex'
set terminal pngcairo size 1600,1000
#set terminal cairolatex pdf color colortext standalone size 20 cm, 6 cm
#set size ratio -1
set xlabel "k/k_F"

set encoding utf8
i=2
rhoBMax=0.09
iRMax=10
rhob=0
set yrange[0:1]

while(i<iRMax){

set output sprintf("form_factor%02d.png", i)
rhoB=rhoBMax-i*rhoBMax/iRMax
set title sprintf("%f fm^{-3}", rhoB)
set multiplot layout 1,2
set key bottom left Left samplen 10

set ylabel offset 0.0,0.0 "F_d^2" rotate by 90

plot  sprintf("../data/form_factor%d.txt", i)  		u 1:2 	w p ps 1 pt 7		lc rgb "#99004d"   	title "droplet" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:3		w l lw 8  dt 1 	lc rgb "#3333cc"   	title "rod, \\parallel" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:4  	w l lw 8  dt 2 	lc rgb "#3333cc"   	title "rod, \\perp" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:5		w l lw 8  dt 1 	lc rgb "#cca300"   	title "slab, \\parallel" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:6		w l lw 8  dt 2 	lc rgb "#cca300"   	title "slab, $\\perp$"

set ylabel offset 0.0,0.0 "\\langle F_d^2\\rangle" rotate by 90
unset title
plot  sprintf("../data/form_factor%d.txt", i)   	u 1:7		w l lw 8  dt 2 		lc rgb "#3333cc"   	title "rod, \\parallel" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:8  	w l lw 15 dt 3 		lc rgb "#3333cc"   	title "rod, \\perp" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:9		w l lw 8  dt 2 		lc rgb "#cca300"   	title "slab, \\parallel" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:10	w l lw 15 dt 3 	lc rgb "#cca300"   	title "slab, \\perp",\
			sprintf("../data/form_factor%d.txt", i)   	u 1:($7+2*$8)  	w l lw 15 dt 1 		lc rgb "#3333cc"   	title "rod, \\perp" ,\
			sprintf("../data/form_factor%d.txt", i)   	u 1:($9+2*$10)	w l lw 8 dt 1 		lc rgb "#cca300"   	title "slab, \\parallel"

i= i+1
unset multiplot
}

