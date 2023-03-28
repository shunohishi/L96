#!/bin/csh
#set compiler="ifort"
set compiler="gfortran"
set module="mod_parameter.f90"
set subroutine="sub_ensemble.f90 sub_initialization.f90 sub_lorenz96.f90 sub_random.f90 sub_obs_operator.f90 sub_sort.f90 sub_lyapunov_equation.f90 sub_etkf.f90 sub_etkfcc.f90 sub_inverce_matrix.f90 sub_eigen.f90 sub_static.f90 sub_iocc.f90 sub_cor_obs_err.f90"
set option="-O3 -llapack"
#set option="-O3 -mkl"
#set debug=""
set debug="-g -fbacktrace"

rm -f lorenz96_enkfcc.out
rm -r ETKFCC
mkdir ETKFCC
${compiler} ${module} lorenz96_enkfcc.f90 ${subroutine} -o lorenz96_enkfcc.out ${option} ${debug}
./lorenz96_enkfcc.out
rm -f *.mod lorenz96_enkfcc.out
