# Gnuplot file
set title 'Harmonic oscillator'
set xlabel 'Time (s)'
set ylabel 'Position (m)'
set grid
set output 'oscil1.eps'
set term postscript eps "Times-Roman" 24
plot 'oscil1.res' u 1:2 t '' w l
