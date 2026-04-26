# Hydrogen ground-state wavefunction comparison.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_wavefunction.svg'

set title 'Hydrogen ground state: numeric vs exact wavefunction'
set xlabel 'r'
set ylabel 'f(r)'
set grid
set key right top

plot 'hydrogen_wavefunction.data' using 1:2 with lines lw 2 title 'numeric', \
     'hydrogen_wavefunction.data' using 1:3 with lines lw 2 dt 2 title 'exact r*exp(-r)'
