set output 'transition.tex'
set terminal cairolatex pdf color colortext standalone size 35 cm, 20 cm

set multiplot layout 2,2 rowsfirst 

set xlabel offset 0.0,0.0 "\\rho_B (fm^{-3})"

# set xtics 0.2

#set xrange[:0.13]
#set yrange[2:10]

set key right top Left box opaque spacing 1.5 font ",10"
set border back

#set rmargin 20
C= "0.680000"
sqrtD= "130.000000"
temp1="0.000000"
parametrization="gm1"
# set yrange[1e-4:1]
# set xrange[0.02:]

set ylabel offset 0.0,0.0 "$\\mu_B$ (MeV)"

plot "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 1:5 w l lw 8 dt 1 lc rgb "black"   title "Hadrons" ,\
      "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 2:6 w l lw 8 dt 2 lc rgb "blue"    title "Quarks" ,\

set ylabel offset 0.0,0.0 "P (MeV fm$^{-3}$)"

# set yrange[0:]

plot "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 1:3 w l lw 8 dt 1 lc rgb "black"   title "Hadrons" ,\
      "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 2:4 w l lw 8 dt 2 lc rgb "blue"    title "Quarks" ,\

set ylabel offset 0.0,0.0 "P (MeV fm$^{-3}$)"
set xlabel offset 0.0,0.0 "$\\mu_B$ (MeV)"
set xrange[900:]
plot "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 5:3 w l lw 8 dt 1 lc rgb "black"   title "Hadrons",\
     "../data/eos_H_".parametrization."_Q_".C."_".sqrtD."_T".temp1.".txt"   	u 6:4 w l lw 8 dt 2 lc rgb "blue"   title "Quarks"

unset multiplot