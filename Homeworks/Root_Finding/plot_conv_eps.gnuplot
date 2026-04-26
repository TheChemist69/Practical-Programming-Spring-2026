# Convergence of E0 with respect to ODE relative tolerance eps.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_conv_eps.svg'

set title 'Hydrogen shooting convergence vs ODE eps'
set xlabel 'eps'
set ylabel '|E0 + 1/2|'
set grid
set logscale xy
set key off

plot 'hydrogen_convergence_eps.data' using 1:3 with linespoints lw 2 pt 7
