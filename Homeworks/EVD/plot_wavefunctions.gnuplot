set terminal svg enhanced size 900,700 background rgb "white"
set output "plot_wavefunctions.svg"
set title "Hydrogen s-wave reduced radial wave functions f_k(r)"
set xlabel "r [Bohr radii]"
set ylabel "f(r)"
set grid
set key top right
set xrange [0:*]

plot \
  "wavefunctions.data" using 1:2 with lines lw 2 title "f_1(r) numerical", \
  "wavefunctions.data" using 1:5 with lines lw 2 dt 2 title "f_1(r) exact", \
  "wavefunctions.data" using 1:3 with lines lw 2 title "f_2(r) numerical", \
  "wavefunctions.data" using 1:6 with lines lw 2 dt 2 title "f_2(r) exact", \
  "wavefunctions.data" using 1:4 with lines lw 2 title "f_3(r) numerical", \
  "wavefunctions.data" using 1:7 with lines lw 2 dt 2 title "f_3(r) exact"
