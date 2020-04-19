set grid
set title 'Simulation of a DC motor - Step response'
set xlabel 'Time (s)'
set ylabel 'Rotation speed (rad/s)' 
plot 'dcmotor.res' u 1:6 w l
