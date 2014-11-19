set terminal png size 800,600
set output 'metis/plot.png'
plot 'metis/speedup.dat' with linespoints
