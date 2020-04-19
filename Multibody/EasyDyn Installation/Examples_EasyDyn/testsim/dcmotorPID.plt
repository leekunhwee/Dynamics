set title 'Simulation of a controlled DC motor'
set xlabel 'Time (s)'
set ylabel 'Motor speed (rad/s)'
set grid
plot 'dcmotorPID.res' u 1:6 w l
