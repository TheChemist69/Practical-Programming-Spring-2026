set terminal svg enhanced size 900,1050 background rgb "white"
set output "plot_cubic.svg"

set multiplot layout 3,1 title "Cubic Spline Interpolation of cos(x)" font ",14"

set grid
set key top right

set title "Cubic interpolant vs cos(x) -- gnuplot csplines comparison"
set xlabel "x"
set ylabel "y"
plot "cubic_interp.txt" using 1:2 with lines lw 2 lc rgb "#0072BD" title "Cubic spline", \
     "cubic_interp.txt" using 1:3 with lines lw 2 dt 2 lc rgb "#D95319" title "cos(x)", \
     "cos_points.txt"   using 1:2 smooth csplines lw 2 dt 4 lc rgb "#77AC30" title "gnuplot csplines", \
     "cos_points.txt"   using 1:2 with points pt 7 ps 1.5 lc rgb "#000000" title "Data points"

set title "Derivative of cubic spline vs -sin(x)"
set xlabel "x"
set ylabel "y'"
plot "cubic_interp.txt" using 1:4 with lines lw 2 lc rgb "#0072BD" title "Spline derivative", \
     "cubic_interp.txt" using 1:5 with lines lw 2 dt 2 lc rgb "#D95319" title "-sin(x)"

set title "Integral of cubic spline vs sin(x)"
set xlabel "x"
set ylabel "Integral"
plot "cubic_interp.txt" using 1:6 with lines lw 2 lc rgb "#0072BD" title "Spline integral", \
     "cubic_interp.txt" using 1:7 with lines lw 2 dt 2 lc rgb "#D95319" title "sin(x)"

unset multiplot
