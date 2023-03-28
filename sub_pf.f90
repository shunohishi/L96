!---------------------------------------------------------------------
! Particle Filter |
!---------------------------------------------------------------------
!
! Obs. pdf: Gaussian
! Deterministic algorithm
!
!---------------------------------------------------------------------

subroutine pf(it,nens,nobs,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,R,H,Neff)

  use parameter
  implicit none

  integer iens,iobs,ix

  integer,intent(in) :: it
  integer,intent(in) :: nens,nobs

  double precision,intent(in) :: xf(nx,nens)                        !Forecast
  double precision,intent(inout) :: xfmean(nx),xfsprd(nx)
  double precision,intent(out) :: xa(nx,nens),xamean(nx),xasprd(nx) !Analysis
  double precision,intent(out) :: Neff            !Effective ensmble size 

  double precision,intent(in) :: obs(nobs)    !Observation
  double precision,intent(in) :: R(nobs,nobs) !Observation error covariance matrix
  double precision,intent(in) :: H(nobs,nx)   !Observation operator
  
  double precision Hxf(nobs,nens)  !Hxf
  double precision inv(nobs,nens)  !Innovation
  double precision invT(nens,nobs)  
  double precision Ri(nobs,nobs)   !R^(-1)
  
  double precision Rinv(nobs)      !R^(-1)(y-Hxf)
  double precision l(nens)         !log Likelihood: -0.5* (y-Hxf)T R^(-1) (y-Hxf)
  double precision phi(nens)       !exp(l(i)-min(l))
  double precision w(nens)         !Weight: phi/sum(phi)
  double precision cw(0:nens)      !Cumlulative weight

  double precision E(nens,nens)    !Transform function
  
  !Hxf
  Hxf(:,:)=matmul(H(:,:),xf(:,:))

  !Innovation
  do iens=1,nens
     inv(:,iens)=obs(:)-Hxf(:,iens)
  end do
  invT(:,:)=transpose(inv(:,:))

  !R^(-1)
  Ri(:,:)=R(:,:)
  call inverce_matrix(nobs,Ri)

  !log p(y/x(i))
  do iens=1,nens
     
     Rinv(:)=matmul(Ri(:,:),inv(:,iens))
     l(iens)=-0.5d0*dot_product(invT(iens,:),Rinv(:))
          
  end do
  
  do iens=1,nens
     phi(iens)=exp(l(iens)-minval(l))
  end do
  
  do iens=1,nens
     w(iens)=phi(iens)/sum(phi)
  end do
  
  !Cumulative sum beta
  cw(0)=0.d0
  do iens=1,nens
     cw(iens)=cw(iens-1)+w(iens)
  end do
  
  if(execute_pf == 1)then !SU
     call stochastic_universal(it,nens,cw,E)
  else if(execute_pf == 2)then !MN
     call multinominal(it,nens,cw,E)
  else if(execute_pf == 3)then !Residual
     call residual_resampling(it,nens,w,E)
  else
     write(*,*) "***Error: Choose execute_pf=1-4"
     stop
  end if
  
  !Assign particle
  xa(:,:)=matmul(xf(:,:),E(:,:))
  
  !Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))     
  end do

  !Effective ensemble size
  Neff=0.d0
  do iens=1,nens
     Neff=Neff+w(iens)**2.d0
  end do
  Neff=1.d0/Neff
  
end subroutine pf

