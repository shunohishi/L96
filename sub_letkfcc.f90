subroutine letkfcc(nens,nobs,sigma_loc,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,Ro,H,inf_parm,cor,inan)

  use parameter
  implicit none

  integer ix,iobs,iens,jens

  integer,intent(in) :: nens,nobs
  integer,intent(in) :: sigma_loc !Localization scale [grid]
  
  !---Global
  double precision,intent(in) ::  xf(nx,nens)           !Ensemble forecast
  double precision,intent(out) :: xfmean(nx),xfsprd(nx) !Ensemble forecast mean/sprd
  double precision,intent(out) :: xa(nx,nens)           !Ensemble analysis
  double precision,intent(out) :: xamean(nx),xasprd(nx) !Ensemble analysis mean/sprd 
  double precision,intent(in) ::  obs(nobs)             !Observation
  
  double precision,intent(in) :: Ro(nobs,nobs)          !Original observation error matrix
  double precision,intent(in) :: H(nobs,nx)             !Observation operator
  double precision,intent(in) :: inf_parm               !Inflation parameter (Multiplicative, RTPP)
  double precision,intent(in) :: cor                    !Correlation coefficient

  integer,intent(inout) :: inan
  
  double precision xpos(nx),opos(nobs)                  !Position

  double precision Xfp(nx,nens),Xap(nx,nens)            !Ensemble perturbation
  
  double precision Yf(nobs,nens)              !yb(iens)=Yf(iens)
  double precision Yfp(nobs,nens)             !Yf=HXf
  double precision Yfmean(nobs),Yfsprd(nobs)  !yfmean = Yfmean

  double precision Yfp_c(nobs,nens)                !Correlated ensemble forecast: (sigma^o_cor/sigma^f(i))*cor*Yfp(i)
  double precision Yfp_oc(nobs,nens)               !Yfp(i)-Yfp_c(i) = [ 1-cor*(sigma^o_cor/sigma^f(i)) ]*Yfp(i)

  double precision r                              ! (sigma^o_cor/sigma^o_orig) * sqrt(1-cor^2)

  !--Adaptive multiplicative inflation
  double precision inv(nobs)
  double precision Roinv(nobs,nobs)
  double precision HPfcH(nobs,nobs)
  double precision rho_u,rho_b
  double precision rho                         !Multiplicative inflation

  
  !---Local
  integer nlobs !Local observation number
  double precision dist,dist_loc
  double precision,allocatable :: lYfp_oc(:,:),lYfmean(:)
  double precision,allocatable :: lobs(:)
  double precision,allocatable :: lRo(:,:),linvRo(:,:)
  double precision,allocatable :: C(:,:)       !C=(Yf)T R-1
  double precision I(nens,nens)                !Unit vector
  double precision Pa(nens,nens)               !Pa
  double precision LWa(nens,nens)              !Wa=sqrt[(nens-1)Pa]
  double precision swa(nens)                   !wa=PaC(y-yb)
  double precision Lswa(nens,nens)             !wa(i)=Wa+wa
  double precision evalue(nens),evector(nens,nens) !Eigen value, vector
  
  double precision inc(nx,nens)                !Increment
  
  !---Position
  do ix=1,nx
     xpos(ix)=dble(ix)
  end do
  opos(:)=matmul(H(:,:),xpos(:))
  
  !---Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
  end do

  if(inan /= 0)then
     xfmean(:)=0.d0
     xfsprd(:)=0.d0
     xa(:,:)=0.d0
     xamean(:)=0.d0
     xasprd(:)=0.d0
     return
  end if
  
  !---1.
  !yb(iens) = Yf(iens)
  do iens=1,nens
     Yf(:,iens)=matmul(H(:,:),xf(:,iens))
  end do

  !---2.
  !ybmean = Yfmean
  do iobs=1,nobs
     call ave_std(nens,Yf(iobs,:),Yfmean(iobs),Yfsprd(iobs))
  end do

  !Yf Perturbation (Yb=HXb)
  do iobs=1,nobs
     Yfp(iobs,:)=Yf(iobs,:)-Yfmean(iobs)
     Yfp_c(iobs,:)=err_obs*cor*Yfp(iobs,:)/Yfsprd(iobs)
     Yfp_oc(iobs,:)=Yfp(iobs,:)-Yfp_c(iobs,:)
  end do
  
  !xf perturbation (Xfp)
  do ix=1,nx
     Xfp(ix,:)=xf(ix,:)-xfmean(ix)
  end do

  !Unit vector
  I(:,:)=0.d0
  do iens=1,nens
     I(iens,iens)=1.d0
  end do

  !Adaptive multiplicative inflation
  if(mult_inf)then
     rho=inf_parm
  else
     rho=1.d0
  end if
  
  !  r=err_obs/err_obs*sqrt(1.d0-cor**2.d0)
  r=sqrt(1.d0-cor**2.d0)
  
  !---3. Choose local from global
  dist_loc=2.d0*sqrt(10.d0/3.d0)*sigma_loc !Gaspari and Cohn (1999)
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
     
     allocate(lYfp_oc(nlobs,nens),lYfmean(nlobs))
     allocate(lobs(nlobs),lRo(nlobs,nlobs),linvRo(nlobs,nlobs))
     allocate(C(nens,nlobs))

     !Gather lobs & Calculate lR
     nlobs=0
     lRo(:,:)=0.d0
     do iobs=1,nobs
        dist=min(abs(opos(iobs)-nx-xpos(ix)),abs(opos(iobs)-xpos(ix)))
        dist=min(abs(opos(iobs)+nx-xpos(ix)),dist)
        if(dist < dist_loc)then
           nlobs=nlobs+1
           lYfp_oc(nlobs,:)=Yfp_oc(iobs,:)
           lYfmean(nlobs)=Yfmean(iobs)
           lobs(nlobs)=obs(iobs)
           lRo(nlobs,nlobs)=Ro(iobs,iobs)*exp(0.5d0*(dist/dble(sigma_loc))**2.d0)
        end if
     end do
     
     !---4.
     !R-1
     linvRo(:,:)=lRo(:,:)
     call inverce_matrix(nlobs,linvRo,inan)
     
     !C=(Yf_oc)T Ro-1
     C(:,:)=matmul(transpose(lYfp_oc(:,:)),linvRo(:,:))
     
     !---5.
     !Pa=[(k-1)I/rho + CYf]-1
     Pa(:,:)=r**2.d0*(dble(nens)-1.d0)/rho*I(:,:)+matmul(C(:,:),lYfp_oc(:,:))

     call inverce_matrix(nens,Pa,inan)
     if(inan /= 0)then
        xfmean(:)=0.d0
        xfsprd(:)=0.d0
        xa(:,:)=0.d0
        xamean(:)=0.d0
        xasprd(:)=0.d0
        deallocate(lYfp_oc,lYfmean)
        deallocate(lobs,lRo,linvRo)
        deallocate(C)
        return
     end if
     
     !---6.
     !sqrt[r^2*(k-1)Pa]
     LWa(:,:)=r**2.d0*(dble(nens)-1.d0)*Pa(:,:)

     !Eigenvalue, vector
     call eigen(nens,LWa,evalue,evector,inan)
     if(inan /= 0)then
        xfmean(:)=0.d0
        xfsprd(:)=0.d0
        xa(:,:)=0.d0
        xamean(:)=0.d0
        xasprd(:)=0.d0
        deallocate(lYfp_oc,lYfmean)
        deallocate(lobs,lRo,linvRo)
        deallocate(C)
        return
     end if
     
     !sqrt(evalue) matrix
     LWa(:,:)=0.d0
     do iens=1,nens
        if(evalue(iens) < 0.d0)then
           write(*,'(a,E10.2)') "Nagative eigen value:", evalue(iens)
           inan=-99
           xfmean(:)=0.d0
           xfsprd(:)=0.d0
           xa(:,:)=0.d0
           xamean(:)=0.d0
           xasprd(:)=0.d0
           deallocate(lYfp_oc,lYfmean)
           deallocate(lobs,lRo,linvRo)
           deallocate(C)
           return
        end if
        LWa(iens,iens)=sqrt(evalue(iens))
     end do
     !sqrt(evalue)matrix * evectorT
     LWa(:,:)=matmul(LWa(:,:),transpose(evector(:,:)))
     !Wa = evector * sqrt(evalue)matrix * evectorT
     LWa(:,:)=matmul(evector(:,:),LWa(:,:))

     if(rtpp_inf)then
        !Winf = alpha*I + (1-alpha) * Wtmp
        do jens=1,nens
           do iens=1,nens
              LWa(iens,jens)=(1.d0-inf_parm)*LWa(iens,jens)
              if(iens == jens) LWa(iens,jens)=LWa(iens,jens)+inf_parm
           end do
        end do
     end if
     
     !---7.wa=PaC(y-yb)
     swa(:)=matmul(C(:,:),(lobs(:)-lYfmean(:)))
     swa(:)=matmul(Pa(:,:),swa(:))
     do iens=1,nens
        Lswa(:,iens)=LWA(:,iens)+swa(:)
     end do
  
     !---8.xa(i) = xfmean + Xbwa(i)
     inc(:,:)=matmul(Xfp(:,:),Lswa(:,:)) !Xbwa(i)
     do iens=1,nens
        xa(ix,iens)=xfmean(ix)+inc(ix,iens)
     end do

     deallocate(lYfp_oc,lYfmean)
     deallocate(lobs,lRo,linvRo)
     deallocate(C)
     
  end do
     
  !---9. Skip
  
  do ix=1,nx
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
  end do

  do iens=1,nens
     do ix=1,nx
        if(xa(ix,iens) /= xa(ix,iens))then
           write(*,*) ix,iens,xa(ix,iens)
           stop
        end if
     end do
  end do
  
end subroutine letkfcc
