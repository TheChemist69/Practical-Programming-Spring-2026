set terminal svg enhanced size 800,600 background rgb "white"
set output "plot_times.svg"
set title "Jacobi diagonalization: scaling with matrix size"
set xlabel "Matrix size N"
set ylabel "Time [seconds]"
set grid
set key top left

f(x) = a * x**3
fit f(x) "out.times.data" using 1:2 via a

plot "out.times.data" using 1:2 with linespoints pt 7 ps 1.2 lw 2 title "Measured time", \
     f(x) with lines lw 2 lc rgb "red" title sprintf("Fit: %.2e × N^3", a)
