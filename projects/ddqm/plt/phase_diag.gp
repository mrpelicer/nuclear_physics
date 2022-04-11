set output 'phase_diag.tex'
set terminal cairolatex pdf color colortext standalone size 35 cm, 20 cm

# set multiplot layout 2,2 rowsfirst 


# set xtics 0.2

#set xrange[:0.13]
#set yrange[2:10]

set key right top Left box opaque spacing 1.5 font ",10"
set border back

#set rmargin 20
C= "0.680000"
sqrtD= "130.000000"
parametrization="gm1"
# set yrange[1e-4:1]
# set xrange[0.02:]
set xlabel offset 0.0,0.0 "$\\mu$ (MeV)"
set ylabel offset 0.0,0.0 "$T$ (MeV)"

plot "../data/diagram_H_".parametrization."_Q_".C."_".sqrtD.".txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "Hadrons" ,\

unset multiplot