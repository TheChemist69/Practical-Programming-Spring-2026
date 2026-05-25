set terminal svg size 900,600 enhanced
set output 'plot_erf.svg'
set terminal svg background "white"

set title 'erf(z): adaptive integration vs references'
set xlabel 'z'
set ylabel 'erf(z)'
set grid

plot 'erf_curve.data' using 1:2 with lines lw 2 title 'erf(z) numerical', \
     'erf_curve.data' using 1:3 with lines dt 2 lw 2 title 'std::erf(z)', \
     'erf_tabulated.data' using 1:2 with points pt 7 ps 1.0 lc rgb 'red' title 'tabulated'
