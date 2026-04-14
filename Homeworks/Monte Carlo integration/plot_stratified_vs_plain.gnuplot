set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_stratified_vs_plain.svg'

set title 'Stratified vs plain Monte Carlo (actual error)'
set xlabel 'N'
set ylabel 'actual error'
set grid
set logscale xy
set key left bottom

plot 'stratified_scaling.data' using 1:2 with linespoints lw 2 pt 7 title 'plain MC', \
     'stratified_scaling.data' using 1:3 with linespoints lw 2 pt 5 title 'stratified MC'
