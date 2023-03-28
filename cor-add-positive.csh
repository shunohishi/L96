#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=KF
set ymax=50
set yint1=10
set yint2=10

#---------------------------------------------------------------------------------
# Analysis RMSE |
#---------------------------------------------------------------------------------

set input=${dir}/ave_rmse_fo.dat
set output=cor-add-positive-rmse.ps
set size=6/6
set range=0/1/1/${ymax}
set label=("(a) KF" "(b) KFCC" "(c) (b)-(a)" "(d) Improvement ratio")
set drange=6.5/3/6/0.25

makecpt -T0/1/0.1 -Crainbow -D > rmse.cpt
makecpt -T0/1/0.1 -Cno_green > rmsecontour1.cpt
makecpt -T0/10/1 -Cno_green > rmsecontour2.cpt

makecpt -T-40/40/10 -Cpolar -D > ir.cpt
makecpt -T-100/100/10 -Cno_green > ircontour1.cpt
makecpt -T-500/500/50 -Cno_green > ircontour2.cpt

makecpt -T-0.5/0.5/0.1 -Cpolar -D > dif.cpt
makecpt -T-5/5/0.1 -Cno_green > difcontour1.cpt
makecpt -T-5/5/0.5 -Cno_green > difcontour2.cpt

gawk '{if($8 != 0. && $3 >= 0.) print $3,$2*10000,$8 > "rmsec.20"}' ${input}
gawk '{if($11 != 0. && $3 >= 0.) print $3,$2*10000,$11 > "rmsenc.20"}' ${input}
gawk '{if($8 != 0. && $11 != 0. && $3 >= 0.) print $3,$2*10000,($11-$8)*100./$11 > "ir.20"}' ${input}
gawk '{if($8 != 0. && $11 != 0. && $3 >= 0.) print $3,$2*10000,$8-$11 > "dif.20"}' ${input}

@ i = 1

foreach var(rmsenc rmsec dif ir)

    #BOX
    if($i == 1 || $i == 2)then
	set BA=a0.5f0.1/a${yint2}f${yint1}:"Inflation\040parameter\040(\040\32710@+-5@+)":
    else
	set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a${yint2}f${yint1}:"Inflation\040parameter\040(\040\32710@+-5@+)":
    endif

    if($i == 1)then
	set BA=${BA}WSne
    else if($i == 2)then
	set BA=${BA}wSne
    else if($i == 3)then
	set BA=${BA}WSne
    else if($i == 4)then
	set BA=${BA}wSne
    endif

    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X2 -Y20 -K -P > ${output}
    else if($i % 2 == 0)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X9.5 -K -O >> ${output}    
    else
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X-9.5 -Y-8 -K -O >> ${output}
    endif

    if(${var} == "rmsenc" || ${var} == "rmsec")then
	set color=rmse
	set contour=rmse
	set dA="+a0+gwhite"
    else
	set color=${var}
	set contour=${var}
	set dA="-"
    endif

    
    #Color
    pscontour ${var}.20 -JX -R -C${color}.cpt -I -K -O >> ${output}

    #Black contour
    pscontour ${var}.20 -JX -R -C${contour}contour1.cpt -W2 -A- -K -O >> ${output}
    pscontour ${var}.20 -JX -R -C${contour}contour2.cpt -W4 -A${dA} -G4c -K -O >> ${output}
    
    pstext -JX -R -N -K -O <<EOF >> ${output}
0 52 14 0 0 LB $label[$i]
EOF

    if(${var} == "rmsenc" || ${var} == "rmsec")then
	set dBA=a0.5f0.1:"RMSE":
	set dE=f0.5c
    else if(${var} == "dif")then
	set dBA=a0.5f0.1:"RMSE":
	set dE=0.5c
    else if(${var} == "ir")then
	set dBA=a20f10:"Improvement\040ratio\040(%)":
	set dE=0.5c
    endif

    psscale -D${drange} -B${dBA} -C${color}.cpt -E${dE} -K -O >> ${output}
    
    @ i++
    
end

rm -f *.20 *.cpt

#----------------------------------------------------------------------------
# Pa |
#-----------------------------------------------------------------------------

set input=${dir}/ave_Pa.dat
set output=cor-add-positive-Pa.ps
set size=6/6
set range=0/1/1/${ymax}
set int1=0.1/1
set int2=0.2/2
set label=("(a) KF" "(b) KFCC" "(c) (b)-(a)")
set drange=6.5/3/6/0.25

makecpt -T0/0.15/0.01 -Crainbow -D > Pa.cpt
makecpt -T0/1/0.01 -Cno_green > Pacontour1.cpt
makecpt -T0/1/0.05 -Cno_green > Pacontour2.cpt

makecpt -T-0.06/0.06/0.01 -Cno_green > dif.cpt
makecpt -T-1/1/0.01 -Cno_green > difcontour1.cpt
makecpt -T-5/5/0.03 -Cno_green > difcontour2.cpt

gawk '{if($4 != 0.) print $3,$2*10000,$4 > "Pa.20"}' ${input}
gawk '{if($6 != 0.) print $3,$2*10000,$6 > "Pac.20"}' ${input}
gawk '{if($8 != 0.) print $3,$2*10000,$8 > "Panc.20"}' ${input}
gawk '{if($6 != 0. && $8 != 0.) print $3,$2*10000,$6-$8 > "dif.20"}' ${input}

@ i = 1

foreach var(Panc Pac dif)

    #BOX
    set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a${yint2}f${yint1}:"Inflation\040parameter\040(\040\32710@+-5@+)":
    if(${var} == "dif")then
	set BA=${BA}WSne
    else
	set BA=${BA}Wsne
    endif

    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y20 -K -P > ${output}
    else
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -Y-7.5 -K -O >> ${output}
    endif

    #Color/Contour
    if(${var} == "dif")then
	set color=dif
	set contour=dif
    else
	set color=Pa
	set contour=Pa
    endif
    pscontour ${var}.20 -JX -R -C${color}.cpt -I -K -O >> ${output}
    pscontour ${var}.20 -JX -R -C${contour}contour1.cpt -W2 -A- -K -O >> ${output}
    pscontour ${var}.20 -JX -R -C${contour}contour2.cpt -W4 -A+a0+gwhite -G6c -K -O >> ${output}

    pstext -JX -R -N -K -O <<EOF >> ${output}
0 52 14 0 0 LB $label[$i]
EOF

    if(${var} == "dif")then
	set dBA=a0.03f0.01:"Pa":
    	psscale -D${drange} -B${dBA} -C${color}.cpt -E0.5 -O >> ${output}
    else
	set dBA=a0.05f0.01:"Pa":
    	psscale -D${drange} -B${dBA} -C${color}.cpt -Ef0.5c -K -O >> ${output}
    endif
	
    @ i++
    
end

rm -f *.20 *.grd *.cpt
