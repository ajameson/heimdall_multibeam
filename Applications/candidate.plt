
#set terminal pngcairo transparent enhanced font "arial,10" size 1280, 960
#set output 'candidate.png'


# Set up 4 subplots
TM = 0.06
BM = 0.14
LM = 0.06
RM = 0.06
#DX = 0.06
#DY = 0.14
NX = 2
NY = 2
SX = 1./NX - (LM + RM)
SY = 1./NY - (TM + BM)
#set size 1.0, 1.0 # Note: This doesn't seem to do anything
set multiplot

set bmargin 0; set tmargin 0; set lmargin 0; set rmargin 0

# SNR vs. DM plot#s (vs. filter)
# -----------------------------
#set xrange[0:1000]
set autoscale x
set xlabel "DM [pc cm^{-3}]"
#set ylabel "SNR (approx)"
set ylabel ""
unset key
set ytics auto
set grid noxtics nomxtics ytics mytics lt 9 lw 0.2

PX = 1
PY = 1
set size SX, SY

# TODO: Add 2 more subplots for filter +1 and -1

set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM
plot "snr_dm.dat" u 2:3 w l



# SNR vs. time plots (vs. filter)
# -------------------------------
PX = 0
PY = 0
dx = 0.0
dy = 0.0
sx = SX
sy = SY/4
set bmargin dy; set tmargin dy; set lmargin dx; set rmargin dx
set size sx, sy

set autoscale x
set yrange [-3:10] # TODO: Needs to be set per candidate?
set xlabel "Time [s]"
#set ylabel "SNR (approx)"
set ylabel ""
unset key
#set xtics 0.002
set ytics 2 format ""
set grid noxtics nomxtics ytics mytics lt 9 lw 0.2

set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM + 0*sy
plot 0 w l lt 9, "snr_time.dat" u 2:6 w l lt 1
set xlabel ""
set xtics format ""
set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM + 1*sy
plot 0 w l lt 9, "snr_time.dat" u 2:4 w l lt 1
set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM + 2*sy
plot 0 w l lt 9, "snr_time.dat" u 2:3 w l lt 1
set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM + 3*sy
plot 0 w l lt 9, "snr_time.dat" u 2:5 w l lt 1



# Freq. vs. time plot
# -------------------
PX = 1
PY = 0
set size SX, SY
set origin PX*(SX+LM+RM) + LM, PY*(SY+TM+BM) + BM

set autoscale x
set autoscale y
unset colorbox
unset xtics
unset ytics
unset x2tics
unset y2tics
set xlabel "Time"
set ylabel "Frequency"

# TODO: Work these out
t0 = -0.40
#t1 = t0 + 0.524288
t1 = t0 + 1.048576
f0 = 1581.
f1 = 1181.
dm = 944

k = 4.148808e3

#f(x) = dm * k * (1/x**2 - 1/f0**2)
#g(x) = f0+f1-sqrt(1. / (x / (dm*k) + 1./f0**2))
g(x) = sqrt(1. / (x / (dm*k) + 1./f0**2))

set x2range [t0:t1]
set y2range [f0:f1]

# TODO: Work out equation for DM curve in pixel coords
plot "freq_time.dat" matrix with image, \
     g(x) w l lc 2 axes x2y2



unset multiplot

pause -1
