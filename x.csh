#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=${argv[1]}
    
#---------------------------------------------------------------
# x |
#---------------------------------------------------------------

echo "x"

set input=${dir}/x.dat
set output=${dir}/x.ps
set size=5/5
set range=1/40/0/720
#if(${dir} == "KF" || ${dir} == "3DVAR") set range=1/40/0/720
#if(${dir} == "ETKF") set range=1/40/0/7200
set int1=1/0.25
set int2=2/0.5
set drange=2.5/-2/5/0.25h

set obs=40
set mult=1.01

gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $5 != 0.) print $3,$4,$5 > "rxt.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $6 != 0.) print $3,$4,$6 > "rx.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $7 != 0.) print $3,$4,$7 > "rxf.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $8 != 0.) print $3,$4,$8 > "rxa.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $5 != 0.) print $3,$4,$5-$5 > "dxt.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $6 != 0.) print $3,$4,$6-$5 > "dx.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $7 != 0.) print $3,$4,$7-$5 > "dxf.20"}' ${input}
gawk -v nobs=${obs} -v nmult=${mult} '{if($1 == nobs && $2 == nmult && $8 != 0.) print $3,$4,$8-$5 > "dxa.20"}' ${input}

@ i = 1

foreach type(r d)

    if(${type} == "r")then
	makecpt -T-10/10/1 -Cno_green -D > color.cpt
	set dBA=a5f1:"x":
    else if(${type} == "d")then
	makecpt -T-1/1/0.2 -Cpolar -D > color.cpt
	set dBA=a1f0.2:"x-x@+true@+":
    endif

    foreach var(xt x xf xa)

	set BA=a5f1:"i":/a60f10:"Time\040(day)":
#	if(${dir} == "KF" || ${dir} == "3DVAR") set BA=a5f1:"i":/a60f10:"Time\040(day)":
#	if(${dir} == "ETKF") set BA=a5f1:"i":/a720f60:"Time\040(day)":
	if(1 <= $i && $i <= 3)then
	    set BA=${BA}Wsne
	else if(5 <= $i && $i <= 7)then
	    set BA=${BA}wsne	
	else if($i == 4)then
	    set BA=${BA}WSne
	else if($i == 8)then
	    set BA=${BA}wSne
	endif

	if($i == 1)then
	    psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y22 -K -P > ${output}
	else if($i == 5)then
	    psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X6 -Y16.5 -K -O >> ${output}
	else
	    psbasemap -JX${size} -R${range} -B${BA} -Gwhite -Y-5.5 -K -O >> ${output}
	endif

	if(${type} == "d" && ${var} == "xt")then

	else

	    surface ${type}${var}.20 -G${type}${var}.grd -R${range} -I${int1}
		
	    psmask ${type}${var}.20 -JX -R -I${int2} -K -O >> ${output}
	    grdview ${type}${var}.grd -JX -R -Ccolor.cpt -Qs -K -O >> ${output}
	    psmask -C -K -O >> ${output}

	endif
	if($i % 4 == 0) psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -K -O >> ${output}
	
	@ i++

    end
end

#gv ${output} &
rm *.20 *.grd color.cpt
