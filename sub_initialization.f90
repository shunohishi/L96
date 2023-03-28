!-----------------------------------------------------------------------------
! Initialization |
!-----------------------------------------------------------------------------

subroutine initialization(amp,deg,noise,x)

  use parameter
  implicit none
  
  integer ix

  double precision,intent(in) :: amp,deg !Amplitude,Degree
  double precision,intent(in) :: noise(nx) !Noise
  double precision,intent(out) :: x(nx)

  !---Sin
  do ix=1,nx
     x(ix)=amp*dsin(dble(nwave)*2.d0*pi*dble(ix)/dble(nx) + deg*pi/180.d0)
  end do

  !Noise
  do ix=1,nx
     x(ix)=x(ix)+noise(ix)
  end do
  
end subroutine initialization

!---------------------------------------------------------------
! Initialization of Error Covariance Matrix |
!---------------------------------------------------------------

subroutine initialization_ECM(nobs,R,Pf)

  use parameter
  implicit none

  integer ix,iobs
  integer,intent(in) :: nobs
  
  double precision,intent(out) :: R(nobs,nobs) !Observation error covarinance matrix: R
  double precision,intent(out) :: Pf(nx,nx)    !Forecast error covariance matrix: P
 
  R(:,:)=0.d0
  Pf(:,:)=0.d0
 
  !---Observation error covariance matrix
  do iobs=1,nobs
     R(iobs,iobs)=err_obs*err_obs
  end do
  
  !---Forecast error covariance matrix
  do ix=1,nx
     Pf(ix,ix)=err_fct*err_fct
  end do
  
end subroutine initialization_ECM
