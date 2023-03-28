!----------------------------------------------------------
! Kalman Filter |
!----------------------------------------------------------
!
! * Observation Operator: H = I
! * Nondiagonal elements = 0 in observation error covariance matrix
!----------------------------------------------------------

subroutine kalman_filter(nobs,xf,xa,Pf,Pa,obs,R,H,HT,C)

  use parameter
  implicit none
  
  integer ix,iobs
  integer,intent(in) :: nobs
  
  double precision,intent(in) ::  xf(nx)    !Forecast:    xf
  double precision,intent(in) ::  obs(nobs)   !Observation: y
  double precision,intent(out) :: xa(nx)    !Analysis:    xa
  double precision inv(nobs)                  !Innovation:  y-Hxf
  
  double precision,intent(inout) ::  Pf(nx,nx) !Forecast error covariance matrix: Pf
  double precision,intent(out) :: Pa(nx,nx)    !Analysis error covariance matrix: Pa
  double precision,intent(in) ::  R(nobs,nobs) !Observation error covarinance matrix: R
    
  double precision PfHT(nx,nobs)          !PfHT
  double precision HPfHT(nobs,nobs)       !HPfHT
  double precision RpHPfHT(nobs,nobs)     !R+HPfHT-HC-(HC)^T --> (R+HPfHT-HC-HC^T)^(-1)
  double precision PfHTmC(nx,nobs)        !PfHT-C
  double precision H(nobs,nx),HT(nx,nobs) !Observation operator
  double precision I(nx,nx) !Unit matrix
  double precision K(nx,nobs) !Kalman gain
  double precision KH(nx,nx)  !KH
  double precision inc(nx)  !Increment
  double precision ImKH(nx,nx) !I-KH

  double precision C(nx,nobs) !forecast and obs. error covariance matrix (ef * eo^T)
  double precision HC(nobs,nobs)
  
  !Check diaonal Pf
  !  do ix=1,nx
  !     if(Pf(ix,ix) <= 0.)then
  !        write(*,*) "Negative Pf:",ix,Pf(ix,ix)
  !     end if
  !  end do
  
  !Unit matrix
  I(:,:)=0.d0
  do ix=1,nx
     I(ix,ix)=1.d0
  end do

  !HC
  HC(:,:)=matmul(H(:,:),C(:,:))
  
  !---Analysis estimate
  !Innovation
  inv(:)=obs(:)-matmul(H(:,:),xf(:)) !y-Hx
    
  !Kalman gain
  PfHT(:,:)=matmul(Pf(:,:),HT(:,:))   !PfHT
  HPfHT(:,:)=matmul(H(:,:),PfHT(:,:)) !HPfHT
  
  RpHPfHT(:,:)=R(:,:)+HPfHT(:,:) &    !R+HPfHT-HC-(HC)T
       & -HC(:,:)-transpose(HC(:,:))    
  call inverce_matrix(nobs,RpHPfHT)   !(R+HPfHT)^(-1)
  PfHTmC(:,:)=PfHT(:,:)-C(:,:)        !PfHT-C
  K(:,:)=matmul(PfHTmC(:,:),RpHPfHT(:,:))   !(PfHT-C)(R+HPfHT-HC-(CH)^T)^(-1)
  
  inc(:)=matmul(K(:,:),inv(:))                !K(y-Hxf)
  
  !Analysis
  xa(:)=xf(:)+inc(:) !xa=xf+K(y-Hxf)
  
  !---Analysis error covariance matrix
  if(execute_da == 1)then !KF

     KH(:,:)=matmul(K(:,:),H(:,:))
     
     ImKH(:,:)=I(:,:)-KH(:,:) !I-KH
     
     !Analysis error covariance matrix
     Pa(:,:)=matmul(ImKH(:,:),Pf(:,:)) & !Pa=(I-KH)Pf+KCT
          & +matmul(K(:,:),transpose(C(:,:)))
          
  else if(execute_da == 2)then !3DVAR
     
     Pa(:,:)=Pf(:,:)
     
  end if
  
end subroutine kalman_filter
