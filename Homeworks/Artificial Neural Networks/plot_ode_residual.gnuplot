set terminal svg enhanced size 800,600 background rgb "white"
set terminal svg enhanced size 900,700 enhanced font 'Helvetica,12'
set output 'plot_ode_residual.svg'
set multiplot layout 2,1 title 'Exercise 3: ODE residual and solution error'
set grid
set key top right
set xlabel 'x'
set ylabel 'Phi(y)'
plot 'ann_ode_residual.data' using 1:2 with lines lw 2 lc rgb '#d95f02' title 'Phi[y]', \
     'ann_ode_residual.data' using 1:3 with lines lw 2 lc rgb '#1b9e77' title 'Phi[y]^2'
set xlabel 'x'
set ylabel 'abs error'
plot 'ann_ode_residual.data' using 1:4 with lines lw 2 lc rgb '#7570b3' title '|Fp(x)-sin(x)|'
unset multiplot
