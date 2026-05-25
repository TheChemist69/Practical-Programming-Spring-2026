set terminal svg size 900,600 enhanced
set output 'plot_erf_acc.svg'
set terminal svg background "white"

set title 'erf(1) absolute error vs acc (eps=0)'
set xlabel 'acc'
set ylabel 'error'
set grid
set logscale xy
set key left bottom

plot 'erf_acc_sweep.data' using 1:(($3>0)?$3:1/0) with linespoints lw 2 pt 7 title '|erf(1)-exact|', \
     'erf_acc_sweep.data' using 1:(($4>0)?$4:1/0) with linespoints lw 2 pt 5 title 'estimated error'
