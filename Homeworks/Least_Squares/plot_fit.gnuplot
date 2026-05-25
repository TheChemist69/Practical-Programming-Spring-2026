set terminal svg enhanced size 900,700 background rgb "white"
set output "plot_fit.svg"
set title "ThX (^{224}Ra) radioactive decay — Rutherford & Soddy (1902)\nOLS fit: y(t) = a·exp(−λt),  T_{1/2} = ln2/λ"
set xlabel "Time [days]"
set ylabel "Activity [relative units]"
set grid
set key top right

plot \
  "fit.txt" using 1:3:4 with filledcurves lc rgb "#ccddff" title "uncertainty band (c ± δc)", \
  "fit.txt" using 1:2 with lines lw 2 lc rgb "#0000cc" title "best fit: a·exp(−λt)", \
  "data.txt" using 1:2:3 with yerrorbars pt 7 ps 1.2 lc rgb "red" title "data"
