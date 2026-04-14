set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_plain_error.svg'

set title 'Plain Monte Carlo: estimated and actual error vs N'
set xlabel 'N'
set ylabel 'error'
set grid
set logscale xy
set key left bottom

f(x) = a/sqrt(x)
a = 1.0
fit f(x) 'plain_scaling.data' using 1:3 via a

plot 'plain_scaling.data' using 1:2 with linespoints lw 2 pt 7 title 'estimated error', \
     'plain_scaling.data' using 1:3 with linespoints lw 2 pt 5 title 'actual error', \
     f(x) with lines lw 2 dt 2 title sprintf('fit a/sqrt(N), a=%.3g', a)
