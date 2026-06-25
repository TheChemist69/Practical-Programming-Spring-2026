set terminal svg size 820,520 background "white"
set output "plot_bessel_convergence.svg"
set title "Contour J_0(2.5): error vs number of polygon sides (Cauchy independence)"
set xlabel "polygon sides"
set ylabel "|J_0 contour - reference|"
set logscale y
set grid
set key top right
plot "bessel_polygon_convergence.data" u 1:2 w lp pt 7 ps 1.0 lw 2 lc rgb "#1f77b4" t "absolute error"
