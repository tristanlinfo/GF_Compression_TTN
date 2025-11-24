# ==============================================================
# Gnuplot script for visualizing TTN cross-interpolation results
# ==============================================================

set terminal pdfcairo enhanced font "Helvetica,14" linewidth 2 size 7,4
set output "ttn_results.pdf"

set style data linespoints
set key outside top center horizontal box
set grid lw 1 lt rgb "#aaaaaa"

set xlabel "Rank"
set ylabel "Mean absolute error"
set logscale y
set format y "10^{%L}"

# --- Define colors for different topologies ---
set style line 1 lc rgb "#1f77b4" lw 2 pt 7  ps 1.2 # BTTN
set style line 2 lc rgb "#ff7f0e" lw 2 pt 9  ps 1.2 # CTTN
set style line 3 lc rgb "#2ca02c" lw 2 pt 5  ps 1.2 # QTT_Block
set style line 4 lc rgb "#d62728" lw 2 pt 11 ps 1.2 # QTT_Alt

# ==============================================================
# 1. Plot error vs bond dimension
# ==============================================================

set title "Mean absolute error vs Bond dimension (G^{<})"
plot \
    "data/bonderrorlist.dat" using 1:2 title "BTTN"       with linespoints ls 1, \
    "" using 1:3 title "CTTN"       with linespoints ls 2, \
    "" using 1:4 title "QTT_Block"  with linespoints ls 3, \
    "" using 1:5 title "QTT_Alt"    with linespoints ls 4

set output "ttn_results_bond.pdf"

replot

# ==============================================================
# 2. Plot error vs number of iterations
# ==============================================================

set output "ttn_results_iter.pdf"
set title "Mean absolute error vs Number of iterations"
set xlabel "Number of iterations"

plot \
    "data/maxitererrorlist.dat" using 1:2 title "BTTN"       with linespoints ls 1, \
    "" using 1:3 title "CTTN"       with linespoints ls 2, \
    "" using 1:4 title "QTT_Block"  with linespoints ls 3, \
    "" using 1:5 title "QTT_Alt"    with linespoints ls 4

set output


# ==============================================================
# Gnuplot script for visualizing TTN cross-interpolation results
# ==============================================================

# --------------------------------------------------------------

set terminal pdfcairo size 10in,6in enhanced font "Helvetica,14"
set output 'results/Chiret_Bond_plot.pdf'
set datafile separator whitespace

set title "Bond-dimension comparison" font ",16"
set xlabel "Bond dimension" font ",14"
set ylabel "Relative error" font ",14"
set key outside right top vertical Right spacing 1.2
set grid lw 0.6 lc rgb "#d0d0d0"
set border lc rgb "#666666"

set style line 1 lt 1 lw 2 pt 7 ps 1.2 lc rgb "#1f77b4" # CTTN (blue)
set style line 2 lt 1 lw 2 pt 9 ps 1.2 lc rgb "#ff7f0e" # QTT_Block (orange)
set style line 3 lt 1 lw 2 pt 11 ps 1.2 lc rgb "#2ca02c" # QTT_Alt (green)

set logscale y
set yrange [1e-10:1]
set format y '10^{%L}'

set xtics nomirror
set ytics nomirror
set mxtics
set mytics

plot 'data/bonderrorlist.dat' using 1:2 with linespoints ls 1 title "CTTN",
'' using 1:3 with linespoints ls 2 title "QTT_Block",
'' using 1:4 with linespoints ls 3 title "QTT_Alt"