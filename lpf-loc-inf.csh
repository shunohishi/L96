#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set input=LPF/ave_bias_rmse.dat
set output=lpf-loc-inf.ps
set size=8/8l
set range=1/10/1e-3/1e0
set BA=a1:"Localization\040scale":/a1:"Additive\040inflation\040factor\040Q":WSne

makecpt -T0.2/1/0.1 -Crainbow -D > color.cpt
makecpt -T0/1/0.1 -Cno_green > contour1.cpt
makecpt -T0/5/1 -Cno_green > contour2.cpt
set drange=9/4/8/0.25
set dBA=a0.2f0.1:"Analysis\040RMSE":

gawk '{if($1 == 512) print $4,$3,$11 > "dat.20"}' ${input}

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}
pscontour dat.20 -JX -R -Ccolor.cpt -I -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour1.cpt -W2 -G8c -A+a0+gwhite -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour2.cpt -W4 -G10c -A+a0+ggray -K -O >> ${output}
psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -K -O >> ${output}

rm dat.20
rm *.cpt
