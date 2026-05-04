set terminal svg enhanced size 900,900 font 'Helvetica,12' background rgb 'white'
set output 'plot_training.svg'

set multiplot layout 2,1 title 'Training cost history for ANN exercises'

# Exercise 1
set title 'Exercise 1: interpolation'
set xlabel 'iteration'
set ylabel 'cost'
set grid
set logscale y
set key top right
plot 'ann_train_history.data' using 1:2 with lines lw 2 lc rgb '#1b9e77' title 'interpolation'


# Exercise 3
set title 'Exercise 3: ODE training'
set xlabel 'iteration'
set ylabel 'cost'
set grid
set logscale y
set key top right
plot 'ann_ode_history.data' using 1:2 with lines lw 2 lc rgb '#d95f02' title 'ODE training'

unset multiplot
