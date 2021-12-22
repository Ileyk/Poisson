set xrange[0:1]
set yrange[0:1]
# set zrange[0:1]
unset autoscale
do for [i=1:100] {p sprintf("output_electrons_%06.0f.dat", i) using 1:3, sprintf("output_positrons_%06.0f.dat", i) using 1:3; pause 0.1 }
