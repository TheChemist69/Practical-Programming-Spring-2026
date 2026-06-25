set terminal svg size 820,520 background "white"
set output "plot_real_trig.svg"
set title "Trig integral by contour:  ∫_0^{2π} dθ/(2 + b cosθ)  vs  2π/√(4 − b²)"
set xlabel "b   (with a = 2 fixed)"
set ylabel "integral value"
set grid
set key top left
plot \
    "real_trig.data" u 1:3 w l lw 3 lc rgb "#bbbbbb" t "exact  2π/√(4 − b²)", \
    "real_trig.data" u 1:2 every 4 w p pt 7 ps 0.8 lc rgb "#1f77b4" t "contour integral"
