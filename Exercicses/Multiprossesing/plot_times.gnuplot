set terminal svg size 800,600 enhanced font 'Arial,14'
set output 'plot.svg'
set title 'Running Time vs Number of Threads'
set xlabel 'Number of threads'
set ylabel 'Time (seconds)'
set grid
set xrange [0:9]
set xtics 1
set key top right
plot 'out.times' using 1:2 with linespoints pointtype 7 pointsize 1.5 lc rgb 'blue' title 'Real time', \
     'out.times' using 1:3 with linespoints pointtype 5 pointsize 1.5 lc rgb 'red' title 'User time'
