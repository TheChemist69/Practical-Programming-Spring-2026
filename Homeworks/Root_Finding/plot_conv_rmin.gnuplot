# Convergence of E0 with respect to rmin.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_conv_rmin.svg'

set title 'Hydrogen shooting convergence vs rmin'
set xlabel 'rmin'
set ylabel '|E0 + 1/2|'
set grid
set logscale xy
set key off

plot 'hydrogen_convergence_rmin.data' using 1:3 with linespoints lw 2 pt 7
