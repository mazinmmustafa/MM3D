reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure1.tex'

set key out vert
set key reverse Left
set key top right
set key spacing 2
# set key box opaque
# set border back

# set xrange [0:200]
# set yrange [-20:0]

# set xtics 10
# set ytics 5

# set mxtics 1
# set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

#set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc #'black' lw 2

set xlabel '$x$ [cm]'
set ylabel '$|E|$ [V/m]'
set title ''

plot 'data1.dat' using 1:2 with lines lw 2 dt 1 lc 'blue' title '$|E_{x}|$', \
     'data1.dat' using 1:3 with lines lw 2 dt 1 lc 'red' title '$|E_{y}|$', \
     'data1.dat' using 1:4 with lines lw 2 dt 1 lc 'black' title '$|E_{z}|$'
