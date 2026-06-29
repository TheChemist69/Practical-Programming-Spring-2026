set terminal svg size 820,520 background "white"
set output "plot_work_precision.svg"
set title "Error estimate and true error vs evaluations  (integral of e^z, 0 to 1+i)" noenhanced
set xlabel "function evaluations (log scale)"
set ylabel "error (log scale)"
set logscale xy
set grid
set key top right noenhanced
eps = 2.22e-16
plot \
    eps w l dt 2 lw 1.5 lc rgb "#888888" t "machine precision (1e-16)", \
    "work_precision.data" u 1:4 w lp pt 5 ps 1.2 lw 2 lc rgb "#d62728" t "estimated error (self-reported)", \
    "work_precision.data" u 1:2 w lp pt 7 ps 1.2 lw 2 lc rgb "#1f77b4" t "true error (vs exact)"
