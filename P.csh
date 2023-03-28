#/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=${argv[1]}

#---------------------------------------------------------------
# P |
#---------------------------------------------------------------

echo "P"

set input=${dir}/P.dat
set output=${dir}/P.ps
set size=6/6
set range=1/40/365/720

set int=1/0.25
set int2=2/0.5

makecpt -T0/0.16/0.02 -Crainbow -D > color.cpt
set drange=7/3/6/0.25
set dBA=a0.5f0.1:"P":

gawk '{if($3 != 0.) print $1,$2,$3 > "Pf.20"}' ${input}
gawk '{if($4 != 0.) print $1,$2,$4 > "Pa.20"}' ${input}

@ i = 1

foreach var(Pf Pa)

    surface ${var}.20 -G${var}.grd -R${range} -I${int1}

    set BA=a5f1:"i":/a60f10:"Time(day)":
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
