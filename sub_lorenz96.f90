!------------------------------------------------------------------------------------------
! Forwading lorenz96 model |
!------------------------------------------------------------------------------------------

subroutine forward_lorenz96(xin,xout)

  use parameter
  implicit none
  
  double precision,intent(in) :: xin(nx)
  double precision,intent(out) :: xout(nx)
  
  double precision xtmp(nx)
  double precision k1(nx),k2(nx),k3(nx),k4(nx)
  
  !---4th-order Runge Kutta
  !k1
  xtmp(:)=xin(:)
  call lorenz96(xtmp,k1)
  
  !k2
  xtmp(:)=xin(:)+0.5d0*k1(:)
  call lorenz96(xtmp,k2)
  
  !k3
  xtmp(:)=xin(:)+0.5d0*k2(:)
  call lorenz96(xtmp,k3)
  
  !k4
  xtmp(:)=xin(:)+k3(:)
  call lorenz96(xtmp,k4)
  
  !xout
  xout(:)=xin(:) + (k1(:) +2.d0*k2(:) +2.d0*k3(:) +k4(:))/6.d0
  
end subroutine forward_lorenz96

!------------------------------------------------------------------------------------------
! lorenz96 model |
!------------------------------------------------------------------------------------------

subroutine lorenz96(xin,xout)

  use parameter
  implicit none
  
  integer ix
  integer ixp1,ixm2,ixm1 !ix+1,ix-2,ix-1

  double precision,intent(in) :: xin(nx)
  double precision,intent(out) :: xout(nx)

  do ix=1,nx
     
     ixp1=ix+1
     ixm2=ix-2
     ixm1=ix-1

     if(nx < ixp1) ixp1=ixp1-nx
     if(ixm2 < 1) ixm2=ixm2+nx
     if(ixm1 < 1) ixm1=ixm1+nx
     
     xout(ix)=(xin(ixp1)-xin(ixm2))*xin(ixm1) -xin(ix) +force
     
  end do

  xout(:)=dt*xout(:)
  
end subroutine lorenz96

!--------------------------------------------------------------------------------------------
