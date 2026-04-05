set terminal svg size 1400,900 background "white"
set output "plot_three_body.svg"
set title "Figure-8 Three-Body Problem"
set xlabel "x"
set ylabel "y"
set size ratio -1
set key outside right top box

# Plot only trajectories. Use a shared dotted pattern and increasing
# linewidth so overlapped sections remain distinguishable.
set style line 1 lc rgb "#2E86AB" lw 1.4 dt (1,3)
set style line 2 lc rgb "#2A9D8F" lw 2.4 dt (1,3)
set style line 3 lc rgb "#7B2CBF" lw 3.4 dt (1,3)

plot \
    "figure8.data" using 2:3 with lines ls 1 title "Body 1 path", \
    "figure8.data" using 4:5 with lines ls 2 title "Body 2 path", \
    "figure8.data" using 6:7 with lines ls 3 title "Body 3 path"
