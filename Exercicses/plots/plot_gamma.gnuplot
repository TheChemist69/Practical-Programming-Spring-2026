set title 'Gamma Function (Stirling Approximation) and Factorials'
set xlabel 'x'
set ylabel 'Gamma(x)'
set grid
set xrange [-5:5]
set yrange [-5:5]
# Add blue dotted vertical lines at x=-4, -3, -2, -1, 0
set style line 10 lt 2 lc rgb 'blue' lw 3 dashtype (2,2)
set arrow from -4, graph 0 to -4, graph 1 nohead ls 10
set arrow from -3, graph 0 to -3, graph 1 nohead ls 10
set arrow from -2, graph 0 to -2, graph 1 nohead ls 10
set arrow from -1, graph 0 to -1, graph 1 nohead ls 10
set arrow from 0, graph 0 to 0, graph 1 nohead ls 10
plot 'gamma_curve.dat' with lines lw 3 title 'Gamma(x) approx', \
     'gamma_points.dat' with points pointtype 7 pointsize 1.5 lw 3 lc rgb 'red' title 'Tabulated factorials'
pause -1 'Press Enter to exit'
