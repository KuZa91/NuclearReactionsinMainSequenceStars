set terminal x11
set autoscale
set title "CNO Reactions"
unset key
plot "NuclearAbundances.txt" using 1:5 with lines title "Y12"
replot "NuclearAbundances.txt" using 1:6 with lines title "Y14"
replot "NuclearAbundances.txt" using 1:7 with lines title "Y16"
