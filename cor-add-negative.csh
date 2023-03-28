#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=KF
set output=cor-add-negative.ps
set size=6/6
set range=-1/0/1/100
set label=("(a) Analysis RMSE" "(b) diagonal-mean Pa")
set drange=3/-2.5/6/0.25h

#----------------------------------------------------------
# Analysis RMSE |
#----------------------------------------------------------

makecpt -T0/0.25/0.025 -Crainbow -D > color.cpt
makecpt -T0/1/0.025 -Cno_green > contour1.cpt
makecpt -T0/1/0.1 -Cno_green > contour2.cpt
set dBA=a0.1f0.025:"Analysis\040RMSE":

set input=${dir}/ave_rmse_fo.dat
gawk '{if($11 != 0. && $3 <= 0.) print $3,$2*10000,$11 > "dat.20"}' ${input}
set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a20f10:"Inflation\040parameter\040(\040\32710@+-5@+)":
set BA=${BA}WSne

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}
pscontour dat.20 -JX -R -Ccolor.cpt -I -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour1.cpt -W2 -A- -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour2.cpt -W4 -A+a0+gwhite  -G2c -K -O >> ${output}
pstext -JX -R -N -K -O <<EOF >> ${output}
-1 105 14 0 0 LB $label[1]
EOF
psscale -D${drange} -B${dBA} -Ccolor.cpt -Ef0.5c -K -O >> ${output}

#---------------------------------------------------------
# Pa |
#---------------------------------------------------------

makecpt -T0/0.3/0.01 -Crainbow -D > color.cpt
makecpt -T0/1/0.01 -Cno_green > contour1.cpt
makecpt -T0/10/0.1 -Cno_green > contour2.cpt
set dBA=a0.1f0.02:"Pa":

set input=${dir}/ave_Pa.dat
gawk '{if($8 != 0.) print $3,$2*10000,$8 > "dat.20"}' ${input}
set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a20f10:"Inflation\040parameter\040(\040\32710@+-5@+)":
set BA=${BA}wSne

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X8 -K -O >> ${output}
pscontour dat.20 -JX -R -Ccolor.cpt -I -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour1.cpt -W2 -A- -K -O >> ${output}
pscontour dat.20 -JX -R -Ccontour2.cpt -W4 -A+a0+gwhite  -G5c -K -O >> ${output}
pstext -JX -R -N -K -O <<EOF >> ${output}
-1 105 14 0 0 LB $label[2]
EOF
psscale -D${drange} -B${dBA} -Ccolor.cpt -Ef0.5c -O >> ${output}

#----------------------------------------------------------------------
# Correlation -0.1 |
#----------------------------------------------------------------------

set output=cor01-add-negative.ps
set size=12/6
set label=("(a) RMSE" "(b) Pa")

set input=${dir}/ave_rmse_fo.dat 
gawk '{if($3 == -0.1) print $2*10000,$8 > "rmsec.20"}' ${input}
gawk '{if($3 == -0.1) print $2*10000,$11 > "rmsenc.20"}' ${input}
set input=${dir}/ave_Pa.dat
gawk '{if($3 == -0.1) print $2*10000,$6 > "Pac.20"}' ${input}
gawk '{if($3 == -0.1) print $2*10000,$8 > "Panc.20"}' ${input}

@ i = 1
foreach var(rmse Pa)

    if(${var} == "rmse")then
	set range=0/100/0.05/0.25
	set BA=a50f10/a0.05f0.01:"RMSE":WSne
    else if(${var} == "Pa")then
	set range=0/100/0/0.15
	set BA=a50f10:"Inflation\040parameter\040(\040\32710@+-5@+)":/a0.05f0.01:"Pa":WSne
    endif

    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}
	pstext -JX -R -N -K -O <<EOF >> ${output}
0 0.27 14 0 0 LB ${label[$i]}
EOF
    else
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -Y-8 -K -O >> ${output}
	pstext -JX -R -N -K -O <<EOF >> ${output}
0 0.16 14 0 0 LB ${label[$i]}
EOF
    endif    

    foreach type(c nc)

	if(${type} == "c") set pen="black"
	if(${type} == "nc") set pen="gray"
	psxy ${var}${type}.20 -JX -R -W5,${pen} -K -O >> ${output}
    
    end

    @ i++
    
end


rm *.20 *.cpt
