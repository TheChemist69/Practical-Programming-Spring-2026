set terminal svg background "white" size 800,600
set output "plot_relativistic.svg"
set title "Relativistic Orbit Precession"
set xlabel "x"
set ylabel "y"
set size ratio -1

plot \
    "orbit_newtonian_circle.data" using (1/$2)*cos($1):(1/$2)*sin($1) with lines title "Newtonian Circular", \
    "orbit_newtonian_ellipse.data" using (1/$2)*cos($1):(1/$2)*sin($1) with lines title "Newtonian Elliptical", \
    "orbit_relativistic.data" using (1/$2)*cos($1):(1/$2)*sin($1) with lines title "Relativistic Precession"
