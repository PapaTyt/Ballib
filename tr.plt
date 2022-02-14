 set linetype  1 lc "blue"    lw 4
 set linetype  2 lc "red"     lw 4
 set linetype  3 lc "magenta" lw 1
 set linetype  4 lc "green"   lw 1
 set linetype  5 lc "gold"    lw 1
 set linetype  6 lc "orange"  lw 1
 
set terminal pngcairo size 1800,1000 enhanced font 'Tahoma,14'
#set xdata time;
#set timefmt "%Y.%m.%d_%H:%M:%S";

#set format x "%d.%m %H:%M";
set grid x y mx my
set key left

#set xtics 21600
set xrange [*:*]

set mxtics 2
set yrange [*:*]
set mytics 2
set xlabel "time" rotate by 0 center offset 0,0,0
set key opaque


f1 = 'tr_full.txt'
f2 = 'tr_M87_full.txt'
f3 = 'tr_SGR_A_full.txt'
set output 'tr_L2_full.png'
set multiplot
r = 2
c = 2
k1=1.0/r
k2=1.0/c

#==================================================
i=0
j=0
set notitle 
#set ylabel "dR[m]" rotate by 0 center offset 0,0,0
set size k2,k1
set origin j*k2,1-(i+1)*k1
#set margin 5,6,3.3
plot f1 using 3:4 w l lt 4 title "X-Y",\
	 f2 using 3:4 w l lt 1 title "X-Y M_8_7",\
	 f3 using 3:4 w l lt 2 title "X-Y SGR-A"
	

#==================================================
i=1
j=0
set notitle 
#set ylabel "dV[mm/s]" rotate by 0 center offset 0,0,0
set size k2,k1
set origin j*k2,1-(i+1)*k1
#set margin 5,6,3.3
plot f1 using 4:5 w l lt 4 title "Y-Z",\
	 f2 using 4:5 w l lt 1 title "Y-Z M_8_7",\
	 f3 using 4:5 w l lt 2 title "X-Y SGR-A"


#==================================================
i=0
j=1
set notitle 
#set ylabel "dR[m]" rotate by 0 center offset 0,0,0
set size k2,k1
set origin j*k2,1-(i+1)*k1
#set margin 5,6,3.3

plot f1 using 3:5 w l lt 4 title "X-Y",\
	 f2 using 3:5 w l lt 1 title "X-Z M_8_7",\
	 f3 using 3:5 w l lt 2 title "X-Y SGR-A"

#==================================================
#i=1
#j=1
#set notitle 
#set ylabel "dV[mm/s]" rotate by 0 center offset 0,0,0
#set size k2,k1
#set origin j*k2,1-(i+1)*k1
#set margin 5,6,3.3

#plot f1 using 3:15 w l lt 15 title "dvx в J2000, мм/с",\
#	 f1 using 3:16 w l lt 16 title "dvy в J2000, мм/с",\
#	 f1 using 3:17 w l lt 17 title "dvz в J2000, мм/с",\
#	 f1 using 3:19 w l lt 2 title "dV в J2000, мм/с"