set terminal svg size 820,520 background "white"
set output "plot_real_rectangle.svg"
set title "∫_{-∞}^{∞} dx/(1+x²) = π:  the non-real edges vanish as the rectangle grows"
set xlabel "rectangle half-width R   (log scale)"
set ylabel "|π − bottom-edge integral|   (log scale)"
set logscale xy
set grid
unset key
plot "real_rectangle.data" u 1:3 w lp pt 7 ps 1.1 lw 2 lc rgb "#d62728"
