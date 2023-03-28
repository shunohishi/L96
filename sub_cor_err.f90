!------------------------------------------------------------------------------------
! Make correlated random number |
!--------------------------------
!
! noise_z=a*noise_x+b*noise_y
!
! noise_x, noise_y: Independent random number with ave = 0.
! sd_x, sd_y, sd_z: standard deviation of noise_x, noise_y, noise_z
! cor_xz: correlation between X and Z
!
!------------------------------------------------------------------------------------

subroutine make_cor_random(n,sd_x,sd_y,sd_z,cor_xz,noise_x,noise_y,noise_z)

  implicit none

  integer i
  integer,intent(in) :: n
  
  double precision,intent(in) :: cor_xz
  double precision,intent(in) :: sd_x,sd_y,sd_z
  double precision,intent(in) :: noise_x(n),noise_y(n)
  double precision,intent(out) :: noise_z(n)

     if(sd_x == 0.d0)then
        noise_z(:)=noise_y(:)
     else if(sd_y == 0.d0)then
        noise_z(:)=noise_x(:)
     else
        do i=1,n
           noise_z(i)=sd_z*(cor_xz*noise_x(i)/sd_x + sqrt(1.d0-cor_xz*cor_xz)*noise_y(i)/sd_y)
        end do
     end if
  
end subroutine make_cor_random

!-----------------------------------------------------------------------------------
! Make correlated observation error |
!------------------------------------------------------------------------------------

subroutine make_cor_err(nens,xt,xf,xferr,obserr,obserr_c,cor)

  use parameter
  implicit none

  !IN
  integer,intent(in) :: nens
  double precision,intent(in) :: xt(nx)
  double precision,intent(in) :: xf(nx,nens)
  double precision,intent(in) :: obserr(nx)
  double precision,intent(in) :: cor

  !OUT
  double precision,intent(out) :: xferr(nx)
  double precision,intent(out) :: obserr_c(nx)

  !TMP
  integer ix,iens
  double precision xfmean(nx),xfsprd(nx)
  double precision xferrmean(nx),xferrsprd(nx)

  !--- Forecast error
  !Ensemble mean/spread
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
  end do

  xferr(:)=xfmean(:)-xt(:)

  if(cor == 0.d0)then
     obserr_c(:)=obserr(:)
     return
  end if
  
  !---Correlated Observation error
  !---Case 1. Forecast spread at each grid 
  if(fo_stdf == 1)then
     
     !Forecast spread from xf(iens)-xt at each grid
     do ix=1,nx

        xferrsprd(ix)=0.d0
        
        do iens=1,nens
           xferrsprd(ix)=xferrsprd(ix)+(xf(ix,iens)-xt(ix))**2.d0
        end do
        
        xferrsprd(ix)=sqrt(xferrsprd(ix)/dble(nens))

     end do
     
  !---Case 2. Spatially averaged forecast spread 
  else if(fo_stdf == 2)then
     
     !Forecast spread from xferr(ix)
     
     call ave_std(nx,xferr,xferrmean(1),xferrsprd(1))
     xferrsprd(:)=xferrsprd(1)
     
  end if

  !Correlated observation error
  do ix=1,nx
     call make_cor_random(nx,xferrsprd(ix),err_obs,err_obs,cor, &
          & xferr(ix),obserr(ix),obserr_c(ix))
  end do

end subroutine make_cor_err
