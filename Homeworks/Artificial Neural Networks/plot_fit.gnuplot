set terminal svg enhanced size 800,600 background rgb "white"
set terminal svg size 900,600 enhanced font 'Helvetica,12'
set output 'plot_fit.svg'
set title 'Exercise 1: ANN interpolation of g(x)'
set xlabel 'x'
set ylabel 'y'
set grid
set key top right
plot 'ann_fit.data' using 1:2 with lines lw 2 lc rgb '#444444' title 'target g(x)', \
     'ann_fit.data' using 1:3 with lines lw 2 lc rgb '#d95f02' title 'ANN F_p(x)'
