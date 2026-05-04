set terminal svg enhanced size 800,600 background rgb "white"
set terminal svg size 900,800 enhanced font 'Helvetica,12'
set output 'plot_derivatives.svg'
set multiplot layout 2,1 title 'Exercise 2: derivatives and antiderivative'
set grid
set key top right
set xlabel 'x'
set ylabel 'derivative'
plot 'ann_derivatives.data' using 1:2 with lines lw 2 lc rgb '#1b9e77' title "target g'(x)", \
     'ann_derivatives.data' using 1:3 with lines lw 2 lc rgb '#d95f02' title "ANN F_p'(x)", \
     'ann_derivatives.data' using 1:4 with lines lw 2 lc rgb '#7570b3' title "target g''(x)", \
     'ann_derivatives.data' using 1:5 with lines lw 2 lc rgb '#e7298a' title "ANN F_p''(x)"
set ylabel 'integral'
set xlabel 'x'
plot 'ann_derivatives.data' using 1:6 with lines lw 2 lc rgb '#1b9e77' title 'target integral', \
     'ann_derivatives.data' using 1:7 with lines lw 2 lc rgb '#d95f02' title 'ANN integral'
unset multiplot
