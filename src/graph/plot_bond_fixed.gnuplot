# Safe plotting script: produce PNG, then PDF (try pdfcairo, fallback to PNG->PDF conversion)
# Usage: gnuplot src/plot_bond_fixed.gnuplot

set datafile separator whitespace

# Common appearance settings
set title "Bond-dimension comparison" font ",16"
set xlabel "Bond dimension" font ",14"
set ylabel "Relative error" font ",14"
set key outside right top vertical Right spacing 1.2
set grid lw 0.6 lc rgb "#d0d0d0"
set border lc rgb "#666666"
set style line 1 lt 1 lw 2 pt 7 ps 1.2 lc rgb "#1f77b4"
set style line 2 lt 1 lw 2 pt 9 ps 1.2 lc rgb "#ff7f0e"
set style line 3 lt 1 lw 2 pt 11 ps 1.2 lc rgb "#2ca02c"
set logscale y
set yrange [1e-16:1]
set format y '10^{%L}'
set xtics nomirror
set ytics nomirror

# Make sure results dir exists (gnuplot can call system())
system("mkdir -p results")

# 1) PNG (always produce a raster fallback)
set terminal pngcairo size 1400,840 enhanced font "Helvetica,14"
set output 'results/GComp_Bond_plot.png'
plot 'data/bonderrorlist.dat' using 1:2 with linespoints ls 1 title "CTTN", \
     '' using 1:3 with linespoints ls 2 title "QTT_Block", \
     '' using 1:4 with linespoints ls 3 title "QTT_Alt"
set output

# 2) Try PDF via pdfcairo (preferred vector output)
set terminal pdfcairo size 10in,6in enhanced font "Helvetica,14"
set output 'results/GComp_Bond_plot.pdf'
# replot uses the previous plot settings
replot
set output

# 3) Fallback: if pdf file is empty or missing, convert PNG -> PDF (macOS sips or imagemagick convert)
# This shell snippet checks file size and converts if needed.
system("bash -lc 'if [ ! -s results/GComp_Bond_plot.pdf ]; then echo "Converting PNG->PDF fallback"; if command -v sips >/dev/null 2>&1; then sips -s format pdf results/GComp_Bond_plot.png --out results/GComp_Bond_plot.pdf; elif command -v convert >/dev/null 2>&1; then convert results/GComp_Bond_plot.png results/GComp_Bond_plot.pdf; else echo "No converter (sips/convert) found; PDF may be empty."; fi; fi'")

# Show summary
print sprintf("Wrote: %s (PNG), %s (PDF)\n", "results/GComp_Bond_plot.png", "results/GComp_Bond_plot.pdf")
