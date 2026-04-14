set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_quasi_vs_plain.svg'

set title 'Quasi-random vs pseudo-random Monte Carlo (actual error)'
set xlabel 'N'
set ylabel 'actual error'
set grid
set logscale xy
set key left bottom

plot 'quasi_scaling.data' using 1:2 with linespoints lw 2 pt 7 title 'plain MC', \
     'quasi_scaling.data' using 1:3 with linespoints lw 2 pt 5 title 'quasi-random MC'
