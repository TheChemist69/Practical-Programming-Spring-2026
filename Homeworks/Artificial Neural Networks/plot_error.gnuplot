set terminal svg enhanced size 800,600 background rgb "white"
set terminal svg enhanced size 900,600 enhanced font 'Helvetica,12'
set output 'plot_error.svg'
set title 'Exercise 1: Interpolation error for ANN'
set xlabel 'x'
set ylabel 'error'
set grid
set key top right
plot 'ann_error.data' using 1:4 with lines lw 2 lc rgb '#d95f02' title 'Fp(x)-g(x)', \
     'ann_error.data' using 1:5 with lines lw 2 lc rgb '#1b9e77' title '|error|'
