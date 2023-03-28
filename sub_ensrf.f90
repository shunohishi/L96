subroutine ensrf(nens,nobs,sigma_loc,xf,xfmean,xfsprd,xa,xamean,xasprd,obs,R,H,mult_parm)

  use parameter
  implicit none

  integer ix,ix1,ix2
  integer iobs,iens

  integer,intent(in) :: nobs,nens
  integer,intent(in) :: sigma_loc !Localization scale [grid]
  
  double precision,intent(in) :: xf(nx,nens)            !Ensemble forecast
  double precision,intent(out) :: xfmean(nx),xfsprd(nx) !Ensemble forecast mean/sprd
  double precision,intent(out) :: xa(nx,nens)           !Ensemble analysis
  double precision,intent(out) :: xamean(nx),xasprd(nx) !Ensemble analysis mean/sprd
  double precision,intent(in) :: obs(nobs)              !Observation
  
  double precision,intent(in) :: R(nobs,nobs)           !Observation error matrix
  double precision,intent(in) :: H(nobs,nx)             !Observation operator
  double precision,intent(in) :: mult_parm              !Multiplicative inflation parameter
  
  double precision invmean(nobs) !Inovation for ensemble mean
  double precision incmean(nx)   !Increment for ensemble mean
  double precision incp(nx,nens) !Increment for ensemble perturbation 
  
  double precision xfp(nx,nens)  !Ensemble forecast perturbation
  double precision Hxfp(nobs,nens)
  double precision xap(nx,nens)  !Ensemble analysis perturbation

  double precision Htmp(1,nx)
  double precision Rtmp(1,1)
  double precision invtmp(1)
  
  double precision Pf(nx,nx)     !Forecast error covariance matrix:Pf
  double precision,allocatable :: PfHT(:,:)    !PfHT
  double precision,allocatable :: HPfHT(:,:)   !HpfHT
  double precision,allocatable :: RpHPfHT(:,:) !R+HPfHT

  double precision,allocatable :: K(:,:)  !Kalman gain
  double precision,allocatable :: Kp(:,:) !Kalman gain for perturbation
    
  double precision,allocatable ::  I(:,:)      !Unit matrix
  double precision alpha         !Kp = alpha * K
  double precision ImKH(nx,nx)   !I-KH

  double precision dYdYTpR(nobs,nobs)
  double precision T(nens,nens)                    !transform matrix: xap = xfp * T
  double precision evalue(nens),evector(nens,nens)

  double precision xamean_tmp(nx)
  double precision xap_tmp(nx,nens)
  double precision xa_tmp(nx,nens)
  
  double precision tmp(nobs,nens)
  
  !---Ensemble mean/sprd
  do ix=1,nx
     call ave_std(nens,xf(ix,:),xfmean(ix),xfsprd(ix))
  end do

  !xfp
  do ix=1,nx
     xfp(ix,:)=xf(ix,:)-xfmean(ix)
  end do

  if(SRF_serial)then
     
     !Allocate
     allocate(PfHT(nx,1),HPfHT(1,1),RpHPfHT(1,1))
     allocate(K(nx,1),Kp(nx,1))
     allocate(I(nx,nx))

     !I
     I(:,:)=0.d0
     do ix=1,nx
        I(ix,ix)=1.d0
     end do

     xa_tmp(:,:)=xf(:,:)
     xamean_tmp(:)=xfmean(:)
     xap_tmp(:,:)=xfp(:,:)

     if(mult_inf)then
        xap_tmp(:,:)=sqrt(mult_parm)*xap_tmp(:,:)
     end if
     
     !---EnSRF
     do iobs=1,nobs
        
        !---Update analytical ensemble mean
        !Inovation
        invmean(:)=obs(:)-matmul(H(:,:),xamean_tmp(:)) !y-xfmean     
        
        !Pf=<xfp,xfpT>
        Pf(:,:)=matmul(xap_tmp(:,:),transpose(xap_tmp(:,:)))/dble(nens-1.d0)
        
        !Pf localization
        if(Pf_loc)then
           call Pf_localization(nx,sigma_loc,Pf)
        end if
        
        !---Update analytical ensemble mean
        Htmp(1,:)=H(iobs,:)
        Rtmp(1,1)=R(iobs,iobs)
        invtmp(1)=invmean(iobs)

        !HPfHT
        PfHT(:,:)=matmul(Pf(:,:),transpose(Htmp(:,:))) !PfHT(nx,1)
        HPfHT(:,:)=matmul(Htmp(:,:),PfHT(:,:))          !HPfHT(1,1)
        RpHPfHT(:,:)=Rtmp(:,:)+HPfHT(:,:)              !R+HPfHT(1,1)
        call inverce_matrix(1,RpHPfHT)                 !(R+HPfHT)^(-1) (1,1)

        K(:,:)=matmul(PfHT(:,:),RpHPfHT(:,:))        !PfHT(R+HPfHT)^(-1) (nx,1)

        incmean(:)=matmul(K(:,:),invtmp(:))

        if(K_loc)then
           call K_localization(nx,sigma_loc,Htmp(1,:),incmean)
        endif
        
        xamean(:)=xamean_tmp(:)+incmean(:)

        !---Update analytical ensemble perturbation
        alpha=1.d0+sqrt(Rtmp(1,1)/(HPfHT(1,1)+Rtmp(1,1)))
        alpha=1.d0/alpha

        Kp(:,:)=alpha*K(:,:)

        ImKH(:,:)=I(:,:)-matmul(Kp(:,:),Htmp(:,:))
        xap(:,:)=matmul(ImKH(:,:),xap_tmp(:,:))     
        
        if(K_loc)then
           incp(:,:)=xap(:,:)-xap_tmp(:,:)
           do iens=1,nens
              call K_localization(nx,sigma_loc,Htmp(1,:),incp(:,iens))
           end do
           xap(:,:)=xap_tmp(:,:)+incp(:,:)
        end if
        
        !---xa=xamean+xap
        do ix=1,nx
           xa(ix,:)=xamean(ix)+xap(ix,:)
        end do

        !---Recalculation (Just in case) 
        do ix=1,nx
           call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
        end do

        do ix=1,nx
           xap(ix,:)=xa(ix,:)-xamean(ix)
        end do


        !Set control variable
        xa_tmp(:,:)=xa(:,:)
        xamean_tmp(:)=xamean(:)
        xap_tmp(:,:)=xap(:,:)
        
     end do !iobs
     
     deallocate(PfHT,HPfHT,RpHPfHT)
     deallocate(K,Kp)
     deallocate(I)

  else

     allocate(PfHT(nx,nobs),HPfHT(nobs,nobs),RpHPfHT(nobs,nobs))
     allocate(K(nx,nobs))
     allocate(I(nens,nens))
     
     I(:,:)=0.d0
     do iens=1,nens
        I(iens,iens)=1.d0
     end do
     
     !---Update analytical ensemble mean
     !inovation
     invmean(:)=obs(:)-matmul(H(:,:),xfmean(:))
     
     !Pf
     Pf(:,:)=matmul(xfp(:,:),transpose(xfp(:,:)))/(nens-1.d0)

     !Multiplicative inflation
     if(mult_inf)then
        Pf(:,:)=mult_parm*Pf(:,:)
     end if

     !Pf localization
     if(Pf_loc)then
        call Pf_localization(nx,sigma_loc,Pf)
     end if
     
     !Kalman gain
     PfHT(:,:)=matmul(Pf(:,:),transpose(H(:,:)))
     HPfHT(:,:)=matmul(H(:,:),PfHT(:,:))
     RpHPfHT(:,:)=R(:,:)+HPfHT(:,:)
     call inverce_matrix(nobs,RpHPfHT)
     K(:,:)=matmul(PfHT(:,:),RpHPfHT(:,:))

     !Increment
     incmean(:)=matmul(K(:,:),invmean(:))
     
     !Analysis
     xamean(:)=xfmean(:)+incmean(:)
     
     !---Update analytical ensemble perturbation
     Hxfp(:,:)=matmul(H(:,:),xfp(:,:))
     dYdYTpR(:,:)=matmul(Hxfp(:,:),transpose(Hxfp(:,:)))+dble(nens-1.d0)*R(:,:) !dYdYT+(nens-1)*R
     call inverce_matrix(nobs,dYdYTpR)     
     tmp(:,:)=matmul(dYdYTpR(:,:),Hxfp(:,:))      !(dYdYT+(nens-1)*R)^-1dY
     T(:,:)=I(:,:)-matmul(transpose(Hxfp(:,:)),tmp(:,:)) !I-dYT(dYdYT+(nens-1)*R)^-1dY
     call eigen(nens,T,evalue,evector)
     T(:,:)=0.d0
     do iens=1,nens
        T(iens,iens)=sqrt(evalue(iens))
     end do
     T(:,:)=matmul(T(:,:),transpose(evector(:,:)))
     T(:,:)=matmul(evector(:,:),T(:,:))

     xap(:,:)=matmul(xfp(:,:),T(:,:))
     
     !---xa=xamean+xap
     do ix=1,nx
        xa(ix,:)=xamean(ix)+xap(ix,:)
     end do

     do ix=1,nx
        call ave_std(nens,xa(ix,:),xamean(ix),xasprd(ix))
     end do
     
     deallocate(PfHT,HPfHT,RpHPfHT)
     deallocate(K)
     deallocate(I)
     
  end if
     
end subroutine ensrf
