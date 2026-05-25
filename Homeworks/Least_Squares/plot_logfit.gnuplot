set terminal svg enhanced size 900,700 background rgb "white"
set output "plot_logfit.svg"
set title "ThX (^{224}Ra) decay in log-space — Rutherford & Soddy (1902)\nOLS fit: ln(y) = ln(a) − λ·t"
set xlabel "Time [days]"
set ylabel "ln(Activity)"
set grid
set key top right

plot \
  "logfit.txt" using 1:3:4 with filledcurves lc rgb "#ccddff" title "uncertainty band (c ± δc)", \
  "logfit.txt" using 1:2 with lines lw 2 lc rgb "#0000cc" title "best fit: ln(a) − λ·t", \
  "logdata.txt" using 1:2:3 with yerrorbars pt 7 ps 1.2 lc rgb "red" title "data"
