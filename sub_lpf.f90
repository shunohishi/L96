subroutine lpf(it,nens,nobs,sigma_loc,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,R,H,Neffmean)

  use parameter
  implicit none

  integer iens,jens,iobs,ix
  
  integer,intent(in) :: it
  integer,intent(in) :: nens,nobs

  double precision,intent(in) :: xf(nx,nens)                        !Forecast
  double precision,intent(inout) :: xfmean(nx),xfsprd(nx)
  double precision,intent(out) :: xa(nx,nens),xamean(nx),xasprd(nx) !Analysis
  double precision,intent(out) :: Neffmean                          !Spatial averaged effective ensmble size
  double precision Neff(nx),Neffsprd
  
  double precision,intent(in) :: obs(nobs)    !Observation
  double precision,intent(in) :: R(nobs,nobs) !Observation error covariance matrix
  double precision,intent(in) :: H(nobs,nx)   !Observation operator

  double precision Hxf(nobs,nens) !Hxf
  double precision inv(nobs,nens) !Innovation
  double precision invT(nens,nobs)
  double precision Ri(nobs,nobs)  !R^(-1)

  double precision Rinv(nobs) !R^(-1)(y-Hxf)
  double precision l(nens)    !log Likelihood: -0.5* invT R^(-1) inv
  double precision phi(nens)  !exp(l(i)-min(l))
  double precision w(nens)    !Weight: phi/sum(phi)
  double precision cw(0:nens) !Cumulative weight

  double precision E(nens,nens) !Transform function

  !Localization
  integer,intent(in) :: sigma_loc !Localization scale [grid]

  integer nlobs
  double precision xpos(nx),opos(nobs)
  double precision dist,dist_loc
  double precision,allocatable :: lobs(:)
  double precision,allocatable :: lHxf(:,:)
  double precision,allocatable :: lR(:,:),lRi(:,:)
  double precision,allocatable :: linv(:,:),linvT(:,:)
  double precision,allocatable :: lRinv(:)
  
  !---Position
  do ix=1,nx
     xpos(ix)=dble(ix)
  end do
  opos(:)=matmul(H(:,:),xpos(:))

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

  
  !---Ensemble spread mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))     
  end do

  !---Main loop for LPF
  dist_loc=2.d0*sqrt(10.d0/3.d0)*sigma_loc
  do ix=1,nx

     !Count nlobs
     nlobs=0
     do iobs=1,nobs
        dist=min(abs(opos(iobs)-nx-xpos(ix)),abs(opos(iobs)-xpos(ix)))
        dist=min(abs(opos(iobs)+nx-xpos(ix)),dist)
        if(dist < dist_loc)then
           nlobs=nlobs+1
        end if
     end do

     if(nlobs == 0) cycle
     
     !Allocate
     allocate(lHxf(nlobs,nens))
     allocate(lobs(nlobs),lR(nlobs,nlobs),lRi(nlobs,nlobs))
     allocate(linv(nlobs,nens),linvT(nens,nlobs))
     allocate(lRinv(nlobs))
     
     !Gather lobs & Calculate lR
     nlobs=0
     lR(:,:)=0.d0
     do iobs=1,nobs
        dist=min(abs(opos(iobs)-nx-xpos(ix)),abs(opos(iobs)-xpos(ix)))
        dist=min(abs(opos(iobs)+nx-xpos(ix)),dist)
        if(dist < dist_loc)then
           nlobs=nlobs+1
           lHxf(nlobs,:)=Hxf(iobs,:)
           lobs(nlobs)=obs(iobs)
           lR(nlobs,nlobs)=R(iobs,iobs)
           linv(nlobs,:)=inv(iobs,:)
           linvT(:,nlobs)=invT(:,iobs)
!           lR(nlobs,nlobs)=R(iobs,iobs)*exp(0.5d0*(dist/dble(sigma_loc))**2.d0)
        end if
     end do
     
     !R^(-1)
     lRi(:,:)=lR(:,:)
     call inverce_matrix(nlobs,lRi)
     
     !log Liklehood
     do iens=1,nens
        lRinv(:)=matmul(lRi(:,:),linv(:,iens))
        l(iens)=-0.5d0*dot_product((linvT(iens,:)),lRinv(:))
     end do
     phi(:)=exp(l(:)-minval(l))
     
     !Weight
     w(:)=phi(:)/sum(phi)

     !Cumulative weight
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

     if(ix == nx+1)then
        do jens=1,nens
           do iens=1,nens
              if(E(iens,jens) /= 0.d0)then
                 write(*,*) iens,jens,E(iens,jens)
              end if
           end do
        end do
     end if
        
     !Assign particle
     xa(ix,:)=matmul(xf(ix,:),E(:,:))
     
     !Effective ensemble size
     Neff(ix)=0.d0
     do iens=1,nens
        Neff(ix)=Neff(ix)+w(iens)**2.d0
     end do
     Neff(ix)=1.d0/Neff(ix)

     !Deallocate
     deallocate(lHxf)
     deallocate(lobs,lR,lRi)
     deallocate(linv,linvT)
     deallocate(lRinv)

     
  end do !ix
  !---End main loop
  
  !Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
  end do
  
  !Spatial mean/sprd Neff
  call ave_std(nx,Neff,Neffmean,Neffsprd)
  
end subroutine lpf
