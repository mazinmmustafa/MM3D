reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure1.tex'

# set key out vert
# set key reverse Left
# set key top right 
# set key spacing 2

set xrange [0:180]
set yrange [-20:+30]

# set xtics 10
# set ytics 0.2

# set mxtics 1
# set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

#set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc #'black' lw 2

set xlabel '$\theta_{s}$ [deg]'
set ylabel '$\sigma_{\theta\theta}$ [dB]'
set title ''

plot 'data1.dat' using 1:2 with lines lw 2 dt 1 lc 'blue' notitle
