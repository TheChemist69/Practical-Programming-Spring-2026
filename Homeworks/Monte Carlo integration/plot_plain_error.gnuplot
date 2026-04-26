# Plot 1: plain Monte Carlo error behavior as sample size increases.
# White background was requested for readable export in reports.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_plain_error.svg'

# Figure labeling.
set title 'Plain Monte Carlo: estimated and actual error vs N'
set xlabel 'N'
set ylabel 'error'
set grid
# Log-log view makes power-law trends easier to inspect.
set logscale xy
set key left bottom

# Reference model: error ~ a / sqrt(N).
f(x) = a/sqrt(x)
a = 1.0
fit f(x) 'plain_scaling.data' using 1:3 via a

# Plot estimated error, actual error, and fitted 1/sqrt(N) curve.
plot 'plain_scaling.data' using 1:2 with linespoints lw 2 pt 7 title 'estimated error', \
     'plain_scaling.data' using 1:3 with linespoints lw 2 pt 5 title 'actual error', \
     f(x) with lines lw 2 dt 2 title sprintf('fit a/sqrt(N), a=%.3g', a)
