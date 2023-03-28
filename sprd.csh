#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=${argv[1]}

#---------------------------------------------------------------
# Sprd
#---------------------------------------------------------------

echo "spread"

set input=${dir}/sprd.dat
set output=${dir}/sprd.ps
set size=6/6
set range=1/40/365/720
#    set range=1/40/365/7200

set int1=1/0.25
set int2=2/0.5

makecpt -T0/0.16/0.02 -Crainbow -D > color.cpt
set drange=7/3/6/0.25
set dBA=a0.5f0.1:"sprd":

gawk '{if($3 != 0.) print $1,$2,$3 > "xfsprd.20"}' ${input}
gawk '{if($4 != 0.) print $1,$2,$4 > "xasprd.20"}' ${input}

@ i = 1

foreach var(xfsprd xasprd)

    surface ${var}.20 -G${var}.grd -R${range} -I${int1}

    set BA=a5f1:"i":/a60f10:"Time(day)":
#    set BA=a5f1:"i":/a720f60:"Time(day)":
    if($i == 1) set BA=${BA}WSne
    if($i == 2) set BA=${BA}wSne
    
    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -K -P > ${output}
    else
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X7 -K -O >> ${output}
    endif

    psmask ${var}.20 -JX -R -I${int2} -K -O >> ${output}
    grdview ${var}.grd -JX -R -Ccolor.cpt -Qs -K -O >> ${output}
    psmask -C -K -O >> ${output}
    if($i == 2) psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -O >> ${output}
    
    @ i++

end

rm *.20 *.grd color.cpt
    

#------------------------

echo "spread mean"

set input=${dir}/sprdmean.dat
set output=${dir}/sprdmean.ps
set size=12/6
#    set range=365/7200/0/0.3
#    set BA=a720f60:"Time(day)":/a0.1f0.02:"Spread":WSne
set range=365/720/0/0.3
set BA=a60f10:"Time(day)":/a0.1f0.02:"Spread":WSne

gawk '{if($3 != 0.) print $1,$2 > "xfsprd.20"}' ${input}
gawk '{if($4 != 0.) print $1,$3 > "xasprd.20"}' ${input}
    
psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}
psxy xfsprd.20 -JX -R -W5,orange -K -O >> ${output}
psxy xasprd.20 -JX -R -W5,cyan -O >> ${output}

rm *.20
    
