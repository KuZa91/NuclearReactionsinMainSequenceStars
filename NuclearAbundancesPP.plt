set terminal x11
set title "PP Reactions"
set autoscale
unset key
plot "NuclearAbundances.txt" using 1:2 with lines title "Y1"
replot "NuclearAbundances.txt" using 1:3 with lines title "Y3"
replot "NuclearAbundances.txt" using 1:4 with lines title "Y4"

