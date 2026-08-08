set term dumb 140 80
set xlabel "time"
plot "outputdir/stats.txt" u 1:2 w l title "uvolavg"
 plot "outputdir/stats.txt" u 1:3 w l title "tau_bot", \
      "outputdir/stats.txt" u 1:4 w l title "tau_top"

