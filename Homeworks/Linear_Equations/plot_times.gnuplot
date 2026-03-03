set terminal svg size 800,600 enhanced font 'Arial,14' background 'white'
set output 'plot.svg'
set title 'QR-decomposition: time vs matrix size N'
set xlabel 'N'
set ylabel 'Time (seconds)'
set xrange [100:1100]
set xtics 100, 100, 1000
set grid

# Fit with c*N^3
f(x) = c * x**3
fit f(x) 'out.times.data' using 1:2 via c

plot 'out.times.data' using 1:2 with points pointtype 7 pointsize 1.5 \
         lc rgb 'blue' title 'Measured time', \
     f(x) with lines lc rgb 'red' lw 2 title sprintf('Fit: %.2e N^3', c)
