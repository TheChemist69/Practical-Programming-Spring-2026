set terminal svg background "white"
set output "plot_oscillator.svg"
set title "Harmonic Oscillator"
set xlabel "x"
set ylabel "y"
plot "oscillator.data" using 1:2 with lines title "y(x)", \
     "oscillator.data" using 1:3 with lines title "y'(x)"
