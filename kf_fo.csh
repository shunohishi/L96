#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

#--------------------------------------------------------------------------
# RMSE
#--------------------------------------------------------------------------

set input=KF/ave_rmse_fo.dat
set output=kf_fo.ps
set size=12/6
set range=0/1/0.0/3.0
set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a0.5f0.1:"Analysis\040RMSE":WSne

gawk '{print $2*0.1,$7 > "rmseac.20"}' ${input}
gawk '{print $2*0.1,$10 > "rmseanc.20"}' ${input}

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y20 -K -P > ${output}

foreach rmse(rmseac rmseanc)

    if(${rmse} == "rmseac") set color="black"
    if(${rmse} == "rmseanc") set color="gray"
    
    psxy ${rmse}.20 -JX -R -Sc0.15 -G${color} -K -O >> ${output}
    psxy ${rmse}.20 -JX -R -W5,${color} -K -O >> ${output}
    
end

pstext -JX -R -N -K -O <<EOF >> ${output}
0 3.2 14 0 0 LB (a) Analysis RMSE
EOF

rm rmse*.20

#------------------------------------------------------------------------
# Correlation |
#-------------------------------------------------------------------------

set size=6/6
set range=0/1/-0.2/1.1
set BA=a0.5f0.1:"Prescribed\040correlation\040coefficient":/a0.5f0.1:"Estimated\040correlation\040coeffieicnt":WSne

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -Y-9 -K -O >> ${output}

#foreach cor(corfnc)
foreach cor(corfc corfnc)

    set input=KF/${cor}_avestd.dat
    gawk '{print $1*0.1,$2,$3 > "diag.20"}' ${input}
    gawk '{print $1*0.1,$4,$5 > "offdiag.20"}' ${input}

    if(${cor} == "corfc") set color="black"
    if(${cor} == "corfnc") set color="gray"

    psxy diag.20 -JX -R -Sc0.15 -G${color} -K -O >> ${output}
    psxy diag.20 -JX -R -Ey/5,${color} -K -O >> ${output}
    psxy diag.20 -JX -R -W5,${color} -K -O >> ${output}
    psxy offdiag.20 -JX -R -Sc0.15 -G${color} -K -O >> ${output}
    psxy offdiag.20 -JX -R -Ey/5,${color} -K -O >> ${output}
    psxy offdiag.20 -JX -R -W5,${color},- -K -O >> ${output}
        
end

pstext -JX -R -N -K -O <<EOF >> ${output}
0 1.2 14 0 0 LB (b) Correlation coefficient
EOF
