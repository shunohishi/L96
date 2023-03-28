#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set EnKF=("LETKF" "EnSRF" "PO")

set output=mult_loc.ps
set size=5/5
set range=1/10/1/1.09

set int1=1/0.01
set int2=2/0.02

makecpt -T0.2/1/0.1 -Crainbow -D > color.cpt
makecpt -T0/1/0.1 -Cno_green > contour1.cpt
makecpt -T0/1/0.5 -Cno_green > contour2.cpt

set drange=5.5/2.5/5/0.25
set dBA=a0.2f0.1:"RMSE":

@ i = 1
foreach DIR(${EnKF})

    set input=${DIR}/ave_bias_rmse.dat
    set BA=a2f1:"Localization\040scale":/a0.02f0.01:"MULT":
    if($i == 1) set BA=${BA}WSne
    if($i == 2) set BA=${BA}wSne
    if($i == 3) set BA=${BA}wSne
    
    gawk '{if($11 != 0.) print $4,$3,$11 > "dat.20"}' ${input}
    surface dat.20 -Gdat.grd -R${range} -I${int1}

    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -X2.5 -Y15 -K -P > ${output}
    else
	psbasemap -JX${size} -R${range} -B${BA} -X5.5 -K -O >> ${output}
    endif
    
    psmask dat.20 -JX -R -I${int2} -K -O >> ${output}
    grdview dat.grd -JX -R -Ccolor.cpt -Qs -K -O >> ${output}
    grdcontour dat.grd -JX -R -Ccontour1.cpt -W2 -K -O >> ${output}
    grdcontour dat.grd -JX -R -Ccontour2.cpt -W4+a0+gwhite -K -O >> ${output}
    psmask -C -K -O >> ${output}

    @ i++
    
end

psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -O >> ${output}

rm dat.20 dat.grd
rm color.cpt contour1.cpt contour2.cpt
