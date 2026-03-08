set terminal svg enhanced size 800,600 background rgb "white"
set output "plot_convergence_dr.svg"
set title "Convergence of ground state energy with grid spacing Δr (r_{max}=10)"
set xlabel "Δr [Bohr radii]"
set ylabel "E_0 [Hartree]"
set grid
set key top left

exact = -0.5

plot "convergence_dr.data" using 1:2 with linespoints pt 7 ps 1.2 lw 2 title "E_0(Δr)", \
     exact with lines lt 2 lw 2 lc rgb "red" title "exact = -0.5 Ha"
