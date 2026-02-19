set title 'Error Function (erf) and Tabulated Values'
set xlabel 'x'
set ylabel 'erf(x)'
set grid
set xrange [-3:3]
set yrange [-1:1]
plot 'erf_curve.dat' with lines title 'erf(x) approx', \
     'erf_points.dat' with points pointtype 7 pointsize 1.5 lc rgb 'red' title 'Tabulated values'
pause -1 'Press Enter to exit'
