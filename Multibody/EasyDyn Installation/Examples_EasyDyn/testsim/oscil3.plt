# Gnuplot file
set title 'Harmonic oscillator simulated with a bad time step'
set xlabel 'Time (s)'
set ylabel 'Position (m)'
set grid
set output 'oscil3.eps'
set term postscript eps "Times-Roman" 24
plot [0:5] 'oscil3.res' u 1:2 t 'simulated' w l,cos(2*pi*x) t 'exact'
