#!/bin/csh
gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set dir=${argv[1]}

#---------------------------------------------------------------
# Bias & RMSE |
#---------------------------------------------------------------

echo "Bias & RMSE"

set input=${dir}/bias_rmse.dat
set output=${dir}/bias_rmse.ps
set size=12/6

gawk '{ print $1,$2 > "xbias.20"}' ${input}
gawk '{ print $1,$3 > "xfbias.20"}' ${input}
gawk '{ print $1,$4 > "xabias.20"}' ${input}
gawk '{ print $1,$5 > "obsbias.20"}' ${input}
gawk '{ print $1,$6 > "xrmse.20"}' ${input}
gawk '{ print $1,$7 > "xfrmse.20"}' ${input}
gawk '{ print $1,$8 > "xarmse.20"}' ${input}
gawk '{ print $1,$9 > "obsrmse.20"}' ${input}

@ i = 1

foreach type(bias rmse)

    if(${type} == "bias" && (${dir} == "KF" || ${dir} == "3DVAR"))then
	set range=365/720/-1/1
	set BA=a60f10:"Time(day)":/a0.5f0.1g99:"Bias":WSne
    else if(${type} == "bias" && (${dir} == "ETKF" || ${dir} == "LETKF"))then
	set range=365/720/-1/1
	set BA=a60f10:"Time(day)":/a0.5f0.1g99:"Bias":WSne
#	set range=365/7200/-1/1
#	set BA=a720f60:"Time(day)":/a0.5f0.1g99:"Bias":WSne
    else if(${type} == "rmse" && (${dir} == "KF" || ${dir} == "3DVAR"))then
	set range=365/720/0/6
	set BA=a60f10:"Time(day)":/a1f0.5:"RMSE":WSne
    else if(${type} == "rmse" && (${dir} == "ETKF" || ${dir} == "LETKF"))then
	set range=365/720/0/6
	set BA=a60f10:"Time(day)":/a1f0.5:"RMSE":WSne
#	set range=365/7200/0/6
#	set BA=a720f60:"Time(day)":/a1f0.5:"RMSE":WSne
    endif

    if($i == 1)then
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y15 -K -P > ${output}
    else
	psbasemap -JX${size} -R${range} -B${BA} -Gwhite -Y-9 -K -O >> ${output}
    endif
    
    foreach var(x xf xa obs)
    
	if(${var} == "x") set color="black"
	if(${var} == "xf") set color="orange"
	if(${var} == "xa") set color="cyan"
	if(${var} == "obs") set color="gray"
	
	psxy ${var}${type}.20 -JX -R -W5,${color} -K -O >> ${output}

    end

    @ i++

end

rm *.20
