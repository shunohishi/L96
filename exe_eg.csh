#!/bin/csh
set compiler="gfortran-mp-9"
set module="mod_parameter.f90"
set subroutine="sub_initialization.f90 sub_lorenz96.f90 sub_random.f90 sub_lyapunov_equation.f90 sub_kalman_filter.f90 sub_inverce_matrix.f90 sub_eigen.f90 sub_static.f90 sub_io.f90"
set option="-O3 -llapack"
set debug=""

rm -f lorenz96_eg.out
${compiler} ${module} lorenz96_eg.f90 -o lorenz96_eg.out ${subroutine} ${option} ${debug}
./lorenz96_eg.out
rm -f *.mod lorenz96_eg.out
