# Plot 4: convergence comparison on the difficult singular benchmark.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_difficult.svg'

# Figure labeling.
set title 'Difficult singular integral: convergence by method'
set xlabel 'N'
set ylabel '|estimate - exact|'
set grid
# Log-log scaling is standard for error-convergence plots.
set logscale xy
set key left bottom

# Columns 5/6/7 contain absolute errors for plain/quasi/stratified.
plot 'difficult_integral.data' using 1:5 with linespoints lw 2 pt 7 title 'plain MC', \
     'difficult_integral.data' using 1:6 with linespoints lw 2 pt 5 title 'quasi-random MC', \
     'difficult_integral.data' using 1:7 with linespoints lw 2 pt 9 title 'stratified MC'
