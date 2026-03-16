set terminal svg enhanced size 900,700 background rgb "white"
set output "plot_linear.svg"

set multiplot layout 2,1 title "Linear Spline Interpolation of cos(x)" font ",14"

set grid
set key top right

set title "Linear interpolant vs cos(x)"
set xlabel "x"
set ylabel "y"
plot "linear_interp.txt" using 1:2 with lines lw 2 lc rgb "#0072BD" title "Linear spline", \
     "linear_interp.txt" using 1:3 with lines lw 2 dt 2 lc rgb "#D95319" title "cos(x)", \
     "cos_points.txt"    using 1:2 with points pt 7 ps 1.5 lc rgb "#000000" title "Data points"

set title "Integral of linear spline vs sin(x)"
set xlabel "x"
set ylabel "Integral"
plot "linear_interp.txt" using 1:4 with lines lw 2 lc rgb "#0072BD" title "Spline integral", \
     "linear_interp.txt" using 1:5 with lines lw 2 dt 2 lc rgb "#D95319" title "sin(x)"

unset multiplot
