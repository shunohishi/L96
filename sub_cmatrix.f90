!------------------------------------------------------------------------------------
! Make C matrix |
!------------------------------------------------------------------------------------

subroutine make_C_matrix_fo(nobs,fo_cc,stdf,stdo,Pf,C)

  use parameter
  implicit none

  integer ix,iobs
  
  integer,intent(in) :: nobs

  double precision,intent(in) :: fo_cc,stdf,stdo
  double precision,intent(in) :: Pf(nx,nx)
  double precision,intent(out) :: C(nx,nobs)
  double precision miss
  
  C(:,:)=0.d0
!  if(fo_cc*stdo < stdf)then
  if(2.d0*fo_cc*stdo < stdf .and. fo_cc*stdo < stdf)then
     C=fo_cc*err_obs/stdf*Pf(:,:)
  end if
  
end subroutine make_C_matrix_fo

!-------------------------------------------

subroutine make_C_matrix_so(nobs,nt,q,eo,C)

  use parameter
  implicit none

  integer ix,iobs,it
  
  integer,intent(in) :: nobs
  integer,intent(in) :: nt

  double precision,intent(in) :: q(nx,nt),eo(nobs,nt)
  double precision,intent(out) :: C(nx,nobs)

  C(:,:)=matmul(q(:,:),transpose(eo(:,:)))
  C(:,:)=C(:,:)/dble(nt)
  
end subroutine make_C_matrix_so

!-----------------------------------------------------------------------------------
! Correlation |
!-----------------------------------------------------------------------------------

subroutine correlation_ef_eo(nt,nobs,ef,eo,cor)

  use parameter
  implicit none

  integer ix,iobs
  integer,intent(in) :: nt,nobs

  double precision,intent(in) :: ef(nt,nx),eo(nt,nobs)
  double precision,intent(out) :: cor(nx,nobs)

  do ix=1,nx
     do iobs=1,nobs
        call correlation(nt,ef(:,ix),eo(:,iobs),cor(ix,iobs))
     end do
  end do
  
end subroutine correlation_ef_eo

!------------------------------------------------------------------------
! Check KH |
!-------------------------------------------------------------------------

subroutine check_KH(fo_cc,stdf,stdo,KH,flag)

  use parameter
  implicit none

  integer ix
  integer,intent(out) :: flag
  
  double precision,intent(in) :: fo_cc
  double precision,intent(in) :: stdf,stdo
  double precision,intent(in) :: KH(nx,nx)
  
  do ix=1,nx
     
     if(KH(ix,ix) < stdf/(stdf-fo_cc*stdo))then
        flag=0
     else
        flag=1
        exit
     end if
     
  end do
  
end subroutine check_KH
