
set terminal svg size 800,600
set terminal svg enhanced size 900,700 background rgb "white"
set output "fit.svg"

set xlabel "Energy [GeV]"
set ylabel "Cross section"
set grid

plot \
  "data.out" using 1:2:3 with yerrorbars title "Data", \
  "fit.out" using 1:2 with lines lw 2 title "Breit-Wigner fit"
