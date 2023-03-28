#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

#---------------------------------------------------------------
# Growth rate |
#---------------------------------------------------------------

set input=GR/growth_rate.dat
set output=GR/growth_rate.ps
set size=6/6
set range=0/20/0/5
set BA=a5f1:"Time(day)":/a1f0.2:"RMSE":WSne
#set range=0/10/0/1
#set BA=a5f1:"Time(day)":/a0.5f0.1:"RMSE":WSne

gawk '{print $1,$4 > "mean.20"}' ${input}
gawk '{if(NR % 10 == 0) print $1,$4,$5 > "std.20"}' ${input}

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}
psxy mean.20 -JX -R -W5,black -K -O >> ${output}
psxy std.20 -JX -R -Ey/black -O >> ${output}

rm mean.20 std.20
