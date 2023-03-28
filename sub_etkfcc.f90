subroutine etkfcc(nens,nobs,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,R,H,pco,inf_parm)

  use parameter
  implicit none

  integer ix,iobs,iens
  integer inan
  
  integer,intent(in) :: nens,nobs
  
  double precision,intent(in) :: xf(nx,nens)            !Ensemble forecast
  double precision,intent(out) :: xfmean(nx),xfsprd(nx) !Ensemble forecast mean/sprd
  double precision,intent(out) :: xa(nx,nens)           !Ensemble analysis
  double precision,intent(out) :: xamean(nx),xasprd(nx) !Ensemble analysis mean/sprd 
  double precision,intent(in) :: obs(nobs)              !Observation

  double precision,intent(in) :: R(nobs,nobs)           !Observation error matrix
  double precision,intent(in) :: H(nobs,nx)    !Observation operator
  double precision,intent(in) :: pco           !Parameter for correlated observation
  double precision,intent(in) :: inf_parm      !Inflation parameter
  
  double precision xfp(nx,nens)                !Ensemble perturbation
  
  double precision Hxf(nobs,nens)              !yb(iens)=Hxf(iens)
  double precision Hxfp(nobs,nens)             !Yf=HXf
  double precision Hxfmean(nobs),Hxfsprd(nobs) !yfmean = Hxfmean

  double precision invR(nobs,nobs)             !R-1
  double precision C(nens,nobs)                !C=(Yf)T R-1

  double precision I(nens,nens)                !Unit vector
  double precision Pa(nens,nens)               !Pa
  double precision inv_Pa(nens,nens)           !Pa^(-1)
  double precision LWa(nens,nens)              !Wa=sqrt[(nens-1)Pa]
  double precision swa(nens)                   !wa=PaC(y-yb)
  double precision Lswa(nens,nens)             !wa(i)=Wa+wa
  double precision evalue(nens),evector(nens,nens) !Eigen value, vector
  double precision rho                         !Multiplicative inflation
  double precision inc(nx,nens)                !Increment

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

  !---3. Skip

  !---4.
  !R-1
  if(diag_err_obs)then
     invR(:,:)=0.d0
     do iobs=1,nobs
        invR(iobs,iobs)=1.d0/R(iobs,iobs)
     end do
  else
     invR(:,:)=R(:,:)
     call inverce_matrix(nobs,invR)
  end if
  
  !C=(Yf)T R-1
  C(:,:)=matmul(transpose(Hxfp(:,:)),invR(:,:))

  !---5.
  I(:,:)=0.d0
  do iens=1,nens
     I(iens,iens)=1.d0
  end do
  
  if(mult_inf)then
     rho=inf_parm
  else
     rho=1.d0
  end if

  ![(k-1)I/rho + CYf]
  Pa(:,:)=I(:,:)*dble(nens-1)/rho+(1.d0-pco)*(1.d0-pco)*matmul(C(:,:),Hxfp(:,:))
  
  !Eigen value/vector
  call eigen(nens,Pa,evalue,evector)
  
  !check evalue
  do iens=1,nens
     if(evalue(iens) < 0.)then
        write(*,*) "***Error: Negative eigenvector ***"
        stop
     end if
  end do

  !Pa=[(k-1)I/rho + CYf]-1
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
  !sqrt(eigenvalue) evectorT
  LWa(:,:)=matmul(LWa(:,:),transpose(evector(:,:)))
  !Wa = evector sqrt(eigenvalue) evectorT
  LWa(:,:)=matmul(evector(:,:),LWa(:,:))
  LWa(:,:)=sqrt(dble(nens)-1.d0)*LWa(:,:)
  
  if(rtpp_inf)then
     !Winf = alpha*I + (1-alpha) * Wtmp
     LWa(:,:)=inf_parm*I(:,:)+(1.d0-inf_parm)*LWa(:,:)
  end if
  
  !---7.wa=PaC(y-yb)
  swa(:)=matmul(C(:,:),(obs(:)-Hxfmean(:)))
  swa(:)=(1.d0-pco)*matmul(Pa(:,:),swa(:))
  do iens=1,nens
     Lswa(:,iens)=LWA(:,iens)+swa(:)
  end do


  !---8.xa(i) = xfmean + Xbwa(i)
  inc(:,:)=matmul(xfp(:,:),Lswa(:,:)) !Xbwa(i)
  do iens=1,nens
     xa(:,iens)=xfmean(:)+inc(:,iens)
  end do

  !---9. Skip
  do ix=1,nx
     call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
  end do
  
end subroutine etkfcc
