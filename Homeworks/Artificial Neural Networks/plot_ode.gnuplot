set terminal svg enhanced size 800,600 background rgb "white"
set terminal svg size 900,700 enhanced font 'Helvetica,12'
set output 'plot_ode.svg'
set multiplot layout 2,1 title 'Exercise 3: ODE training for y prime pime prime prime + y = 0'
set grid
set key top right
set xlabel 'x'
set ylabel 'y'
plot 'ann_ode.data' using 1:2 with lines lw 2 lc rgb '#444444' title 'exact sin(x)', \
     'ann_ode.data' using 1:3 with lines lw 2 lc rgb '#d95f02' title 'ANN solution'
set xlabel 'x'
set ylabel "y'"
plot 'ann_ode.data' using 1:5 with lines lw 2 lc rgb '#7570b3' title "exact cos(x)", \
     'ann_ode.data' using 1:4 with lines lw 2 lc rgb '#1b9e77' title "ANN derivative"
unset multiplot
