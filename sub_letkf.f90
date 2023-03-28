subroutine letkf(nens,nobs,sigma_loc,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,R,H,mult_parm)

  use parameter
  implicit none

  integer ix,iobs,iens

  integer,intent(in) :: nens,nobs
  integer,intent(in) :: sigma_loc !Localization scale [grid]
  
  !---Global
  double precision,intent(in) :: xf(nx,nens)            !Ensemble forecast
  double precision,intent(out) :: xfmean(nx),xfsprd(nx) !Ensemble forecast mean/sprd
  double precision,intent(out) :: xa(nx,nens)           !Ensemble analysis
  double precision,intent(out) :: xamean(nx),xasprd(nx) !Ensemble analysis mean/sprd 
  double precision,intent(in) :: obs(nobs)              !Observation
  
  double precision,intent(in) :: R(nobs,nobs)           !Observation error matrix
  double precision,intent(in) :: H(nobs,nx)             !Observation operator
  double precision,intent(in) :: mult_parm              !Multiplicative inflation
  
  double precision xpos(nx),opos(nobs)                  !Position
  double precision xfp(nx,nens),xap(nx,nens)            !Ensemble perturbation
  
  double precision Hxf(nobs,nens)              !yb(iens)=Hxf(iens)
  double precision Hxfp(nobs,nens)             !Yf=HXf
  double precision Hxfmean(nobs),Hxfsprd(nobs) !yfmean = Hxfmean

  !---Local
  integer nlobs !Local observation number
  double precision dist,dist_loc
  double precision,allocatable :: lHxf(:,:),lHxfp(:,:),lHxfmean(:)
  double precision,allocatable :: lobs(:)
  double precision,allocatable :: lR(:,:),linvR(:,:)
  double precision,allocatable :: C(:,:)       !C=(Yf)T R-1
  double precision I(nens,nens)                !Unit vector
  double precision Pa(nens,nens)               !Pa
  double precision LWa(nens,nens)              !Wa=sqrt[(nens-1)Pa]
  double precision swa(nens)                   !wa=PaC(y-yb)
  double precision Lswa(nens,nens)             !wa(i)=Wa+wa
  double precision,allocatable :: evalue(:),evector(:,:) !Eigenvalue/vector
  
  double precision rho                         !Multiplicative inflation
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
  
  !---1.
  !yb(iens) = Hxf(iens)
  do iens=1,nens
     Hxf(:,iens)=matmul(H(:,:),xf(:,iens))
  end do

  !---2.
  !ybmean = Hxfmean
  do iobs=1,nobs
     call ave_std(nens,Hxf(iobs,:),Hxfmean(iobs),Hxfsprd(iobs))
  end do

  !Hxf Perturbation (Yb=HXb)
  do iobs=1,nobs
     Hxfp(iobs,:)=Hxf(iobs,:)-Hxfmean(iobs)
  end do

  !xf perturbation (xfp)
  do ix=1,nx
     xfp(ix,:)=xf(ix,:)-xfmean(ix)
  end do

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
     
     allocate(lHxf(nlobs,nens),lHxfp(nlobs,nens),lHxfmean(nlobs))
     allocate(lobs(nlobs),lR(nlobs,nlobs),linvR(nlobs,nlobs))
     allocate(C(nens,nlobs))

     !Gather lobs & Calculate lR
     nlobs=0
     lR(:,:)=0.d0
     do iobs=1,nobs
        dist=min(abs(opos(iobs)-nx-xpos(ix)),abs(opos(iobs)-xpos(ix)))
        dist=min(abs(opos(iobs)+nx-xpos(ix)),dist)
        if(dist < dist_loc)then
           nlobs=nlobs+1
           lHxf(nlobs,:)=Hxf(iobs,:)
           lHxfp(nlobs,:)=Hxfp(iobs,:)
           lHxfmean(nlobs)=Hxfmean(iobs)
           lobs(nlobs)=obs(iobs)
           lR(nlobs,nlobs)=R(iobs,iobs)*exp(0.5d0*(dist/dble(sigma_loc))**2.d0)
        end if
     end do
     
     !---4.
     !R-1
     if(diag_err_obs)then
        linvR(:,:)=0.d0
        do iobs=1,nlobs
           linvR(iobs,iobs)=1.d0/lR(iobs,iobs)
        end do
     else
        allocate(evalue(nlobs),evector(nlobs,nlobs))
        call eigen(nlobs,lR,evalue,evector)
        linvR(:,:)=0.d0
        do iobs=1,nlobs
           linvR(iobs,iobs)=1.d0/evalue(iobs)
        end do
        linvR(:,:)=matmul(linvR(:,:),evector(:,:))
        linvR(:,:)=matmul(evector(:,:),linvR(:,:))
        deallocate(evalue,evector)
     endif
     
     !C=(Yf)T R-1
     C(:,:)=matmul(transpose(lHxfp(:,:)),linvR(:,:))
     
     !---5.
     I(:,:)=0.d0
     do iens=1,nens
        I(iens,iens)=1.d0
     end do
     if(mult_inf)then
        rho=mult_parm
     else
        rho=1.d0
     end if

     ![(k-1)I/rho + CYf]
     Pa(:,:)=(dble(nens)-1.d0)/rho*I(:,:)+matmul(C(:,:),lHxfp(:,:))

     !Pa=[(k-1)I/rho + CYf]-1
     allocate(evalue(nens),evector(nens,nens))
     call eigen(nens,Pa,evalue,evector)
     do iens=1,nens
        if(evalue(iens) < 0.)then
           write(*,*) "*** Error: Negative eigen value """
           stop
        end if
     end do

     Pa(:,:)=0.d0
     do iens=1,nens
        Pa(iens,iens)=1.d0/evalue(iens) !^-1
     end do
     Pa(:,:)=matmul(Pa(:,:),transpose(evector(:,:)))
     Pa(:,:)=matmul(evector(:,:),Pa(:,:))
     
     !---6.
     !LWa=sqrt[(k-1)Pa]
     LWa(:,:)=0.d0
     do iens=1,nens
        LWa(iens,iens)=1.d0/sqrt(evalue(iens)) !^-0.5
     end do
     LWa(:,:)=matmul(LWa(:,:),transpose(evector(:,:)))
     LWa(:,:)=matmul(evector(:,:),LWa(:,:))
     LWa(:,:)=sqrt(dble(nens-1))*LWa(:,:)
     
     deallocate(evalue,evector)
     
     !---7.wa=PaC(y-yb)
     swa(:)=matmul(C(:,:),(lobs(:)-lHxfmean(:)))
     swa(:)=matmul(Pa(:,:),swa(:))
     do iens=1,nens
        Lswa(:,iens)=LWA(:,iens)+swa(:)
     end do
  
     !---8.xa(i) = xfmean + Xbwa(i)
     inc(:,:)=matmul(xfp(:,:),Lswa(:,:)) !Xbwa(i)
     do iens=1,nens
        xa(ix,iens)=xfmean(ix)+inc(ix,iens)
     end do

     deallocate(lHxf,lHxfp,lHxfmean)
     deallocate(lobs,lR,linvR)
     deallocate(C)
     
  end do
     
  !---9. Skip
  
  do ix=1,nx
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
  end do
  
end subroutine letkf
