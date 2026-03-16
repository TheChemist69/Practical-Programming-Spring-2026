set terminal svg enhanced size 900,1050 background rgb "white"
set output "plot_quadratic.svg"

set multiplot layout 3,1 title "Quadratic Spline Interpolation of sin(x)" font ",14"

set grid
set key top right

set title "Quadratic interpolant vs sin(x)"
set xlabel "x"
set ylabel "y"
plot "quad_interp.txt" using 1:2 with lines lw 2 lc rgb "#0072BD" title "Quadratic spline", \
     "quad_interp.txt" using 1:3 with lines lw 2 dt 2 lc rgb "#D95319" title "sin(x)", \
     "sin_points.txt"  using 1:2 with points pt 7 ps 1.5 lc rgb "#000000" title "Data points"

set title "Derivative of quadratic spline vs cos(x)"
set xlabel "x"
set ylabel "y'"
plot "quad_interp.txt" using 1:4 with lines lw 2 lc rgb "#0072BD" title "Spline derivative", \
     "quad_interp.txt" using 1:5 with lines lw 2 dt 2 lc rgb "#D95319" title "cos(x)"

set title "Integral of quadratic spline vs 1-cos(x)"
set xlabel "x"
set ylabel "Integral"
plot "quad_interp.txt" using 1:6 with lines lw 2 lc rgb "#0072BD" title "Spline integral", \
     "quad_interp.txt" using 1:7 with lines lw 2 dt 2 lc rgb "#D95319" title "1-cos(x)"

unset multiplot
