set terminal png size 1280, 720 enhanced font "Times-New-Roman,18"
if(!exists("filename")) filename='WOA_10parkinsons150421_0.txt'
if(!exists("outputname")) outputname='ConvWOA.png'
if(!exists("titlename")) titlename='Whale Optimization Algorithm'
set output outputname

set key on
set key bottom right
set auto
set grid
set size 1,1
set yrange[0:1]
set ytics 0, 0.1, 1
#stats filename u 1:4 i 0 name "a"
stats filename name "a"
max_col = a_columns
set xrange[0:a_max_x]

set title titlename
if(max_col<10){
    set  multiplot layout 2,1
        #unset xtics
        unset xlabel
        set ylabel 'Fitness\_fold-1'
        plot filename u 1:4 i 0 w l lw 4 title "Best Fitness" , "" u 1:7 i 0 w l lw 4 title "Worst Fitness"
        unset ylabel
        unset title
        set xlabel 'Generaciones'
        set ylabel 'Fitness\_fold-5'
        plot filename u 1:4 i 4 w l lw 4 title "Best Fitness", "" u 1:7 i 4 w l lw 4 title "Worst Fitness"
    unset multiplot
}else{
    set  multiplot layout 2,1
        #unset xtics
        unset xlabel
        set ylabel 'Fitness\_fold-1'
        plot filename u 1:4 i 0 w l lw 4 title "Best Fitness" , "" u 1:7 i 0 w l lw 4 title "Worst Fitness", "" u 1:10 i 0 w l lw 4 title "allTimeBest"
        unset ylabel
        unset title
        set xlabel 'Generaciones'
        set ylabel 'Fitness\_fold-5'
        plot filename u 1:4 i 4 w l lw 4 title "Best Fitness", "" u 1:7 i 4 w l lw 4 title "Worst Fitness", "" u 1:10 i 4 w l lw 4 title "allTimeBest"
    unset multiplot
}
