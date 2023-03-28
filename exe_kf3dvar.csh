#!/bin/csh
#set compiler="ifort"
set compiler="gfortran"
set module="mod_parameter.f90"
set subroutine="sub_initialization.f90 sub_lorenz96.f90 sub_random.f90 sub_sort.f90 sub_obs_operator.f90 sub_lyapunov_equation.f90 sub_kalman_filter.f90 sub_cmatrix.f90 sub_inverce_matrix.f90 sub_eigen.f90 sub_static.f90 sub_io.f90 sub_check_nan.f90"
#set option="-O3 -mkl"
set option="-O3 -llapack"
set debug=""
#set debug="-g -fbacktrace"

rm -f lorenz96_kf3dvar.out
${compiler} ${module} lorenz96_kf3dvar.f90 -o lorenz96_kf3dvar.out ${subroutine} ${option} ${debug}
./lorenz96_kf3dvar.out
rm -f *.mod lorenz96_kf3dvar.out
