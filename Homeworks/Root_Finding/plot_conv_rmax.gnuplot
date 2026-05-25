# Convergence of E0 with respect to rmax.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_conv_rmax.svg'

set title 'Hydrogen shooting convergence vs rmax'
set xlabel 'rmax'
set ylabel '|E0 + 1/2|'
set grid
set logscale y
set key off

plot 'hydrogen_convergence_rmax.data' using 1:3 with linespoints lw 2 pt 7
