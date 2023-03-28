#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set output=add_pf.ps
set size=5/5
set range=0/29/100/5000

#------------------------------------------------------------------------
# Analysis RMSE |
#------------------------------------------------------------------------

set label=("(a) Stochastic Universal" "(b) Multinominal" "(c) Residual")

makecpt -T0.4/1/0.1 -Crainbow -D > color.cpt
makecpt -T0/1/0.1 -Cno_green > contour1.cpt
makecpt -T0/1/1 -Cno_green > contour2.cpt

set drange=2.5/-2/5/0.25h
set dBA=a0.1:"Analysis\040RMSE":

@ i = 1
foreach type(SU MN R)

    set input=PF_${type}/ave_bias_rmse.dat
    gawk '{print $3*100,$1,$11 > "dat.20"}' ${input}
    set BA=a5f1:"Q\040(\32710@+-3@+)":/a500f100W:"Particle\040size\040":

    if($i == 1)then
	set BA=${BA}WSne
    	psbasemap -JX${size} -R${range} -B${BA} -X3 -Y20 -K -P > ${output}
    else
	set BA=${BA}wSne
	psbasemap -JX${size} -R${range} -B${BA} -X5.5 -K -O >> ${output}
    endif
    pscontour dat.20 -JX -R -Ccolor.cpt -I -K -O >> ${output}
    pscontour dat.20 -JX -R -Ccontour1.cpt -W2 -A- -K -O >> ${output}
    pscontour dat.20 -JX -R -Ccontour2.cpt -W4 -A+a0+gwhite -G5c -K -O >> ${output}

    pstext -JX -R -N -K -O <<EOF >> ${output}
0 5500 14 0 0 LB $label[$i]
EOF
    
    if($i == 2)then
	psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -K -O >> ${output}
    endif
    
    @ i++
    
end

#-----------------------------------------------------------------------
# Effective ensmble size |
#-----------------------------------------------------------------------

set label=("(d) Stochastic Universal" "(e) Multinominal" "(f) Residual")

makecpt -T10/150/10 -Crainbow -D > color.cpt
makecpt -T0/500/10 -Cno_green > contour1.cpt
makecpt -T0/500/50 -Cno_green > contour2.cpt

set drange=2.5/-2/5/0.25h
set dBA=a50f10:"Effective\040ensemble\040size":

@ i = 1
foreach type(SU MN R)

    set input=PF_${type}/Neff.dat
    gawk '{if($3 != 0.) print $3*100,$1,$4 > "dat.20"}' ${input}
    set BA=a5f1:"Q\040(\32710@+-3@+)":/a500f100W:"Particle\040size\040":

    if($i == 1)then
	set BA=${BA}WSne
    	psbasemap -JX${size} -R${range} -B${BA} -X-11 -Y-10 -K -O >> ${output}
    else
	set BA=${BA}wSne
	psbasemap -JX${size} -R${range} -B${BA} -X5.5 -K -O >> ${output}
    endif
    pscontour dat.20 -JX -R -Ccolor.cpt -I -K -O >> ${output}
    pscontour dat.20 -JX -R -Ccontour1.cpt -W2 -A- -K -O >> ${output}
    pscontour dat.20 -JX -R -Ccontour2.cpt -W4 -A+a0+gwhite -G5c -K -O >> ${output}

    pstext -JX -R -N -K -O <<EOF >> ${output}
0 5500 14 0 0 LB $label[$i]
EOF
    
    if($i == 2)then
	psscale -D${drange} -B${dBA} -Ccolor.cpt -E0.5c -K -O >> ${output}
    endif
    
    @ i++
    
end

rm dat.20
rm *.cpt

