!------------------------------------------------------------------
! Make parameter "a"
!------------------------------------------------------------------
!
! eo_cor=a*ef+eo_orig
!
!    eo_cor(orig): correlated (original) observation error
!    ef: forecast error
!
!------------------------------------------------------------------

subroutine make_a(ia,a)

  implicit none
  integer,intent(in) :: ia
  double precision,intent(out) :: a

  a=-1+0.1d0*dble(ia-1)
  
end subroutine make_a

!------------------------------------------------------------------
! Make correlated R |
!------------------------------------------------------------------

subroutine make_Rcor(nobs,iR,R)

  implicit none

  integer iobs
  
  integer,intent(in) :: nobs
  integer,intent(in) :: iR
  double precision,intent(out) :: R(nobs,nobs)

  R(:,:)=0.d0

  do iobs=1,nobs
     R(iobs,iobs)=1.d0+0.1d0*dble(iR-1)
  end do
  
end subroutine make_Rcor
