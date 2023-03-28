#!/bin/cs
gmt set FONT_LABEL=14p,Helvetica,black

set size=7/7
set range=-0.9/0.9/0/1.0
set BAx=a0.3f0.1+l"Correlation coefficient"
set BAy=a0.5f0.1+l"Inflation parameter"
set BAl=WSne
set int1=0.1/0.1
set int2=0.2/0.2
set drange=7.5/0.5+w5/0.25+e0.5
set dBA=a0.5f0.1+l"RMSE"

set input=LETKF/ave_rmse.dat
gawk '{if($4 == 1 && $11  != 0.) print $5,$3,$11 > "dat.20"}' ${input}
#gawk '{if($4 == 1 && $14  != 0.) print $5,$3,$14 > "dat.20"}' ${input}
#gmt surface dat.20 -Gdat.grd -R${range} -I${int1}


gmt makecpt -T0/1/0.1 -Crainbow -D > color.cpt

exit

gmt begin letkfcc_rmse png

    gmt basemap -JX${size} -R${range} -Bx${BAx} -By${BAy} -B${BAl} -X3 -Y20

    gmt psmask dat.20 -I${int2}
    gmt pscontour dat.20 -Ccolor.cpt -I
    #    gmt grdview dat.grd -Ccolor.cpt -Qs
    gmt psmask -C
    
    gmt colorbar -Dx${drange} -Bx${dBA} -Ccolor.cpt
    
gmt end

#rm *.20 *.grd *.cpt
