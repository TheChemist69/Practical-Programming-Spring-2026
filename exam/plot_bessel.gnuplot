set terminal svg size 820,520 background "white"
set output "plot_bessel.svg"
set title "Bessel functions J_n(x) by complex contour integration vs series reference"
set xlabel "x"
set ylabel "J_n(x)"
set grid
set key top right
set samples 400
plot \
    "bessel_contour.data" u 1:3 w l lw 3 lc rgb "#bbbbbb" t "J_0 ref", \
    "bessel_contour.data" u 1:2 every 8 w p pt 7 ps 0.6 lc rgb "#1f77b4" t "J_0 contour", \
    "bessel_contour.data" u 1:5 w l lw 3 lc rgb "#dddddd" t "J_1 ref", \
    "bessel_contour.data" u 1:4 every 8 w p pt 7 ps 0.6 lc rgb "#d62728" t "J_1 contour", \
    "bessel_contour.data" u 1:7 w l lw 3 lc rgb "#eeeeee" t "J_2 ref", \
    "bessel_contour.data" u 1:6 every 8 w p pt 7 ps 0.6 lc rgb "#2ca02c" t "J_2 contour"
