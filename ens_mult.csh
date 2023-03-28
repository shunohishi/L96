gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=$argv[1]

set input=${dir}/ave_bias_rmse.dat
set output=${dir}/ave_bias_rmse.ps
set size=12/12
set range=10/50/1/1.1
set BA=a10f1:"Ensemble\040member":/a0.05f0.01:"Inflation\040factor":WSne

set int1=10/0.01
set int2=20/0.02

set drange=13/6/12/0.25
set dBA=a0.2f0.1:"RMSE":
makecpt -T0.2/1/0.1 -Crainbow -D > color.cpt
makecpt -T0/1/0.1 -Cno_green > contour1.cpt
makecpt -T0/1/0.5 -Cno_green > contour2.cpt

gawk '{if($9 != 0.) print $1,$3,$10 > "dat.20"}' ${input}
surface dat.20 -Gdat.grd -R${range} -I${int1}

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}

psmask dat.20 -JX -R -I${int2} -K -O >> ${output}
grdview dat.grd -JX -R -Ccolor.cpt -Qs -K -O >> ${output}
grdcontour dat.grd -JX -R -Ccontour1.cpt -W2 -A+a0+gwhite -G6c -K -O >> ${output}
grdcontour dat.grd -JX -R -Ccontour2.cpt -W4 -A- -K -O >> ${output}
psmask -C -K -O >> ${output}

psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -O >> ${output}

#rm -f dat.20 dat.grd color.cpt contour1.cpt contour2.cpt
