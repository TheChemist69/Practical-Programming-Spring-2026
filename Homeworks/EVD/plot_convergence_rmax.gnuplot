set terminal svg enhanced size 800,600 background rgb "white"
set output "plot_convergence_rmax.svg"
set title "Convergence of E_0 with r_{max}  (fixed {/Symbol D}r = 0.1 a_0)"
set xlabel "r_{max} [Bohr radii]"
set ylabel "E_0 [Hartree]"
set grid
set key top right

exact = -0.5

plot "convergence_rmax.data" using 1:2 with linespoints pt 7 ps 1.2 lw 2 title "E_0(r_{max})", \
     exact with lines lt 2 lw 2 lc rgb "red" title "exact = -0.5 Ha"
