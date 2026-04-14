set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_difficult.svg'

set title 'Difficult singular integral: convergence by method'
set xlabel 'N'
set ylabel '|estimate - exact|'
set grid
set logscale xy
set key left bottom

plot 'difficult_integral.data' using 1:5 with linespoints lw 2 pt 7 title 'plain MC', \
     'difficult_integral.data' using 1:6 with linespoints lw 2 pt 5 title 'quasi-random MC', \
     'difficult_integral.data' using 1:7 with linespoints lw 2 pt 9 title 'stratified MC'
