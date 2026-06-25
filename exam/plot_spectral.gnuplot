set terminal svg size 820,520 background "white"
set output "plot_spectral.svg"
set title "Periodic trapezoid vs adaptive open-4 on a circular contour (J_3(5))"
set xlabel "function evaluations (log scale)"
set ylabel "absolute error"
set logscale xy
set grid
set key top right
plot \
    "spectral_trapezoid.data" u 1:2 w lp pt 7 ps 1.0 lw 2 lc rgb "#1f77b4" t "periodic trapezoid", \
    "spectral_adaptive.data"  u 1:2 w lp pt 5 ps 1.0 lw 2 lc rgb "#d62728" t "adaptive open-4"
