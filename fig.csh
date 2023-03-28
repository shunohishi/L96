#!/bin/csh

set dir=${argv[1]}

csh x.csh ${dir} && convert -density 200 ${dir}/x.ps ${dir}/x.png
csh bias_rmse.csh ${dir} && convert -density 200 ${dir}/bias_rmse.ps ${dir}/bias_rmse.png

if(${dir} == "KF" || ${dir} == "3DVAR")then
    csh P.csh ${dir} && convert -density 200 ${dir}/P.ps ${dir}/P.png
endif

if(${dir} == "ETKF" || ${dir} == "LETKF")then
    csh sprd.csh ${dir} && convert -density 200 ${dir}/sprd.ps ${dir}/sprd.png && convert -density 200 ${dir}/sprdmean.ps ${dir}/sprdmean.png
endif

wait
