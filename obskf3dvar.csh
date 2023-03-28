gmtset PAPER_MEDIA a4+ LABEL_FONT_SIZE 14

set output=obskf3dvar.ps
set size=12/6
set range=1/40/0/5
set BA=a10f1:"Number\040of\040obs.":/a1f0.2:"RMSE":WSne

psbasemap -JX${size} -R${range} -B${BA} -Gwhite -X3 -Y10 -K -P > ${output}

foreach dir(KF 3DVAR)

    set input=${dir}/ave_bias_rmse.dat
    if($dir == "KF")then
	gawk '{if($2 == 1.01 && $9 != 0.) print $1,$9 > "dat.20"}' ${input}
	set color="red"
    else if($dir == "3DVAR")then
	gawk '{if($2 == 1. && $9 != 0.) print $1,$9 > "dat.20"}' ${input}
	set color="blue"
    endif
    
    psxy dat.20 -JX -R -W5,${color} -K -O >> ${output}

end

rm dat.20
