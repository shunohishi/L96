#!/bin/csh
#set compiler="ifort"
set compiler="gfortran"
set module="mod_parameter.f90"
set subroutine="sub_ensemble.f90 sub_initialization.f90 sub_lorenz96.f90 sub_random.f90 sub_obs_operator.f90 sub_sort.f90 sub_lyapunov_equation.f90 sub_etkf.f90 sub_letkf.f90 sub_letkfcc.f90 sub_ensrf.f90 sub_poenkf.f90 sub_poenkfcc.f90 sub_pf.f90 sub_lpf.f90 sub_localization.f90 sub_inverce_matrix.f90 sub_eigen.f90 sub_resampling.f90 sub_static.f90 sub_io.f90"
set option="-O3 -llapack"
#set option="-O3 -mkl"
#set debug=""
set debug="-g -fbacktrace"

rm -f lorenz96_enkf.out
${compiler} ${module} lorenz96_enkf.f90 ${subroutine} -o lorenz96_enkf.out ${option} ${debug}
./lorenz96_enkf.out
rm -f *.mod lorenz96_enkf.out

