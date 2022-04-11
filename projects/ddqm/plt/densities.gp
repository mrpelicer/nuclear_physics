set output 'densities.tex'
set terminal cairolatex pdf color colortext standalone size 35 cm, 20 cm

set multiplot layout 3,2 rowsfirst 

set xlabel offset 0.0,0.0 "\\rho_B (fm^{-3})"
set ylabel offset 0.0,0.0 "\\rho_i (fm^{-3})"

set xtics 0.2

#set xrange[:0.13]
#set yrange[2:10]

set key Left box opaque at 0.5,1.8 spacing 1.5 font ",10"
set border back

#set rmargin 20
C= "0.000000"
sqrtD= "165.000000"
temp1="0.000000"
temp2="50.000000"
set yrange[1e-4:1]
set xrange[:1]
set log y
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "U, T=0" ,\
  	 "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	u 1:2 w l lw 8 dt 2 lc rgb "black"   title "U, T=50" ,\
     "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	u 1:3 w l lw 8 dt 1 lc rgb "blue"    title "D, T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	u 1:3 w l lw 8 dt 2 lc rgb "blue"    title "D, T=50" ,\
     "../data/eos_".C."_".sqrtD."_T".temp1.".txt"    	u 1:4 w l lw 8 dt 1 lc rgb "red"     title "S, T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	u 1:4 w l lw 8 dt 2 lc rgb "red"     title "S, T=50" ,\
     "../data/eos_".C."_".sqrtD."_T".temp1.".txt"    	u 1:5 w l lw 6 dt 1 lc rgb "purple"   title "e, T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	u 1:5 w l lw 6 dt 2 lc rgb "purple"   title "e, T=50" ,\
     "../data/eos_".C."_".sqrtD."_T".temp1.".txt"    	u 1:6 w l lw 6 dt 1 lc rgb "pink"     title "m, T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	u 1:6 w l lw 6 dt 2 lc rgb "pink"     title "m, T=50"

unset yrange
unset log y
set arrow from 0.0,0.0 to 1.6,0.0 dt 2
set ylabel offset 0.0,00 "P (MeV/fm^{-3})"
set key Left box opaque at 0.4,350 spacing 1.5 font ",12"
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	     u 1:8 w l lw 8 dt 1 lc rgb "black"   title " T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	     u 1:8 w l lw 8 dt 2 lc rgb "black"   title "T=50"
      
unset arrow
unset key
set ylabel offset 0.0,00 "\\mu_B (MeV)"
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	     u 1:9 w l lw 8 dt 1 lc rgb "black"   title " T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	     u 1:9 w l lw 8 dt 2 lc rgb "black"   title "T=50"

set ylabel offset 0.0,00 "\\varepsilon (MeV/fm^{-3})"
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	     u 1:7 w l lw 8 dt 1 lc rgb "black"   title " T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"            u 1:7 w l lw 8 dt 2 lc rgb "black"   title "T=50"

set ylabel offset 0.0,00 "\\varepsilon_{2F}/\\rho_B (MeV)"
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"   	     u 1:11 w l lw 8 dt 1 lc rgb "black"   title " T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"   	     u 1:11 w l lw 8 dt 2 lc rgb "black"   title "T=50"

set ylabel offset 0.0,00 "${\\mathcal F}_{3F}/\\rho_B$ (MeV)"
plot "../data/eos_".C."_".sqrtD."_T".temp1.".txt"             u 1:10 w l lw 8 dt 1      lc rgb "black"   title " T=0" ,\
     "../data/eos_".C."_".sqrtD."_T".temp2.".txt"            u 1:10 w l lw 8 dt 2      lc rgb "black"   title "T=50"  ,\
     "../data/stability_ms100.000000_T".temp1.".txt"              u 7:5 w p pt 6 ps 1.    lc rgb"black"   notitle  ,\
     "../data/stability_ms100.000000_T".temp2.".txt"             u 7:5 w p pt 6 ps 1.    lc rgb"black"   notitle  
unset multiplot