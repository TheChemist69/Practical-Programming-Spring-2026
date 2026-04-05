set terminal svg background "white"
set output "plot_lotka_volterra.svg"
set title "Lotka-Volterra System"
set xlabel "t"
set ylabel "Population"
plot "lotka_volterra.data" using 1:2 with lines title "Prey", \
     "lotka_volterra.data" using 1:3 with lines title "Predator"
