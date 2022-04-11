set output 'pressure.tex'
set terminal cairolatex pdf color colortext standalone size 40 cm, 10 cm

set multiplot layout 1,4 rowsfirst 

set xlabel offset 0.0,0.0 "\\mu_B (fm^{-3})"
set ylabel offset 0.0,0.0 "$P$ (Mev  fm$^{-3}$)"

# set xtics 0.2

#set xrange[:0.13]
#set yrange[2:10]


set key Left box opaque top left spacing 1.5 font ",10"
set border back

#set rmargin 20

set xrange[939:1650]
# set yrange[0:]
par1="l3wr"
par2="fsu2h"
plot "../data/press_H_".par1."_T0.000000.txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "L3$\\omega\\rho$" ,\
    "../data/press_H_".par2."_T0.000000.txt"    	u 1:2 w l lw 8 dt 1 lc rgb "red"   title "FSU2H",\
    "../data/press_Q_0.680000_130.000000_T0.000000.txt"   	u 1:2 w l lw 8 dt 2 lc rgb "black"   title "$C=0.2, \\sqrt(D) = 195$ MeV^2$"

plot "../data/press_H_".par1."_T50.000000.txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "L3$\\omega\\rho$" ,\
    "../data/press_H_".par2."_T50.000000.txt"     u 1:2 w l lw 8 dt 1 lc rgb "red"   title "FSU2H"

plot "../data/press_H_".par1."_T100.000000.txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "L3$\\omega\\rho$" ,\
    "../data/press_H_".par2."_T100.000000.txt"    u 1:2 w l lw 8 dt 1 lc rgb "red"   title "FSU2H"

plot "../data/press_H_".par1."_T150.000000.txt"   	u 1:2 w l lw 8 dt 1 lc rgb "black"   title "L3$\\omega\\rho$" ,\
    "../data/press_H_".par2."_T150.000000.txt"      u 1:2 w l lw 8 dt 1 lc rgb "red"   title "FSU2H"

