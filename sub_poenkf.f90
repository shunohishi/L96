subroutine poenkf(nens,nobs,sigma_loc,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,obserr,R,H,mult_parm)

  use parameter
  implicit none

  integer ix,iobs,iens

  integer,intent(in) :: nobs,nens
  integer,intent(in) :: sigma_loc !Localization scale [grid]
  
  double precision,intent(inout) :: xf(nx,nens)            !Ensemble forecast 
  double precision,intent(out) :: xfmean(nx),xfsprd(nx) !Ensemble forecast meam/sprd
  double precision xfp(nx,nens)                         !Ensemble forecast perturbation
  double precision yfp(nobs,nens)
  
  double precision,intent(out) :: xa(nx,nens)           !Ensemble analysis
  double precision,intent(out) :: xamean(nx),xasprd(nx) !Ensemble analysis mean/sprd
  double precision xap(nx,nens)                         !Ensemble analysis perturbation
  
  double precision,intent(in) :: obs(nobs)        !Observation
  double precision,intent(in) :: obserr(nx,nens)  !Observation ensemble error
  double precision obsens(nobs,nens)              !Ensmble observation
  
  double precision,intent(in) :: R(nobs,nobs)     !Obs. error covariance matrix
  double precision,intent(in) :: H(nobs,nx)       !Obs. Operator
  double precision,intent(in) :: mult_parm        !Multiplicative inflation

  double precision Pf(nx,nx)          !Forecast error covariance matrix
  double precision PfHT(nx,nobs)      !PfHT
  double precision HPfHT(nobs,nobs)   !HPfHT
  double precision RpHPfHT(nobs,nobs) !R+HPfHT

  double precision K(nx,nobs)         !Kalman gain
  double precision inv(nobs)          !Innovation
  double precision inc(nx,nens)       !Increment
  
  !Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
  end do

  !Ensemble forecast perturbation
  do ix=1,nx
     xfp(ix,:)=xf(ix,:)-xfmean(ix)
  end do
  
  !Multiplicative inflation
  if(mult_inf)then
     xfp(:,:)=sqrt(mult_parm)*xfp(:,:)
     do ix=1,nx
        xf(ix,:)=xfmean(ix)+xfp(ix,:)
     end do
!     Pf(:,:)=mult_parm*Pf(:,:)
  end if

  yfp(:,:)=matmul(H(:,:),xfp(:,:))
  
  !Pf
  !Pf(:,:)=matmul(xfp(:,:),transpose(xfp(:,:)))/dble(nens-1.d0)
  
  !Pf localization
  if(Pf_loc)then
     call Pf_localization(nx,sigma_loc,Pf)
  end if
  
  !Kalman gain
  !PfHT(:,:)=matmul(Pf(:,:),transpose(H(:,:)))
  !HPfHT(:,:)=matmul(H(:,:),PfHT(:,:))
  PfHT(:,:)=matmul(xfp(:,:),transpose(yfp(:,:)))
  HPfHT(:,:)=matmul(yfp(:,:),transpose(yfp(:,:)))
  RpHPfHT(:,:)=dble(nens-1)*R(:,:)+HPfHT(:,:)
  call inverce_matrix(nobs,RpHPfHT)

!  write(*,*) RpHPfHT(1,1)
  
  K(:,:)=matmul(PfHT(:,:),RpHPfHT(:,:))

  !Ensmble obs
  obsens(:,:)=matmul(H(:,:),obserr(:,:))
  do iobs=1,nobs
     do iens=1,nens
        obsens(iobs,iens)=obsens(iobs,iens)+obs(iobs)
     end do
  end do
  
  !Ensmble analysis
  do iens=1,nens
     inv(:)=obsens(:,iens)-matmul(H(:,:),xf(:,iens))
     inc(:,iens)=matmul(K(:,:),inv(:))
  end do
  xa(:,:)=xf(:,:)+inc(:,:)

  !Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
  end do
  
!  write(*,*) obs(1)
!  write(*,*) xfmean(1),xfsprd(1)
!  write(*,*) xamean(1),xasprd(1)
  !stop
  
end subroutine poenkf

