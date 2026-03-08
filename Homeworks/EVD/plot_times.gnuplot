set terminal svg enhanced size 800,600 background rgb "white"
set output "plot_times.svg"
set title "Jacobi diagonalization: O(N^3) scaling  (cyclic sweeps, random symmetric matrix)"
set xlabel "Matrix size N"
set ylabel "Time [seconds]"
set grid
set key top left

f(x) = a * x**3
g(x) = b * x**3
fit f(x) "out.times.data" using 1:2 via a
fit g(x) "out.times.data" using 1:3 via b

plot "out.times.data" using 1:2 with linespoints pt 7 ps 1.2 lw 2 title "Standard", \
     "out.times.data" using 1:3 with linespoints pt 5 ps 1.2 lw 2 title "Optimized (upper-triangle)", \
     f(x) with lines lw 2 lc rgb "red" dt 2 title sprintf("Fit: %.2e N^3", a), \
     g(x) with lines lw 2 lc rgb "blue" dt 2 title sprintf("Fit: %.2e N^3", b)
