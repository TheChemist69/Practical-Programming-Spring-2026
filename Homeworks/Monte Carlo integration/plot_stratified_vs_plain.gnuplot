# Plot 3: compare actual error for stratified vs plain Monte Carlo.
set terminal svg size 900,600 enhanced background rgb 'white'
set output 'plot_stratified_vs_plain.svg'

# Figure labeling.
set title 'Stratified vs plain Monte Carlo (actual error)'
set xlabel 'N'
set ylabel 'actual error'
set grid
# Log-log coordinates reveal relative rate and constant differences.
set logscale xy
set key left bottom

# Column 2 = plain actual error, column 3 = stratified actual error.
plot 'stratified_scaling.data' using 1:2 with linespoints lw 2 pt 7 title 'plain MC', \
     'stratified_scaling.data' using 1:3 with linespoints lw 2 pt 5 title 'stratified MC'
