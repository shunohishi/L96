
!-------------------------------------------------------------------
! Lyapunov equation |
!-------------------------------------------------------------------
!
! Calculate Lyapnov equation:Pf(t)=MPf(t-1)MT
!
!--------------------------------------------------------------------

subroutine lyapunov_equation(x,Pf)

  use parameter
  implicit none

  integer,parameter :: itype=1
  integer ix
  integer ixm2,ixm1,ixp1

  double precision,intent(in) :: x(nx)
  double precision,intent(inout) :: Pf(nx,nx) !Forecast error covariance matrix
  
  double precision M(nx,nx)  !Tangent linear model/State transition matrix
  double precision xpdx(nx) !x+dx
  double precision Mx(nx),Mxpdx(nx) !M(x),M(x+dx)
  double precision MT(nx,nx) !MT

  if(execute_da == 1)then !KF
  
     !Tangent linear model(state transition matrix): M
     M(:,:)=0.d0
     do ix=1,nx
        
        if(itype == 1)then
           ixm2=ix-2
           ixm1=ix-1
           ixp1=ix+1
           
           if(ixm2 <= 0) ixm2=nx+ixm2
           if(ixm1 <= 0) ixm1=nx+ixm1
           if(ixp1 > nx) ixp1=ixp1-nx

           M(ix,ixm2)=-1.d0*dt*x(ixm1)
           M(ix,ixm1)=dt*(x(ixp1)-x(ixm2))
           M(ix,ix)=1.d0-dt
           M(ix,ixp1)=dt*x(ixm1)

        else
           
           xpdx(:)=x(:)
           xpdx(ix)=xpdx(ix)+dx_fct
           call forward_lorenz96(x,Mx)
           call forward_lorenz96(xpdx,Mxpdx)
           M(:,ix)=(Mxpdx(:)-Mx(:))/dx_fct

        end if
           
     end do
     MT(:,:)=transpose(M(:,:))
     
     !Lyapunov equation
     Pf(:,:)=matmul(Pf(:,:),MT(:,:)) !PfMT
     Pf(:,:)=matmul(M(:,:),Pf(:,:)) !MPfMT

  else if(execute_da == 2)then !3DVAR

     Pf(:,:)=Pf(:,:)

  end if
     
end subroutine lyapunov_equation

!--------------------------------------------------------------------
! Covariance inflation |
!--------------------------------------------------------------------

subroutine covariance_inflation(Pf,inf_parm)

  use parameter
  implicit none

  integer ix
  double precision,intent(inout) :: Pf(nx,nx)
  double precision,intent(in) :: inf_parm
  
  if(execute_da == 1)then !KF
  
     !Multiplicative inflation
     if(mult_inf)then
        !Pf(:,:)=inf_parm*Pf(:,:)
        do ix=1,nx
           Pf(ix,ix)=inf_parm*Pf(ix,ix)
        end do
     end if
     
     !Additive inflation
     if(add_inf)then
        do ix=1,nx
           Pf(ix,ix)=inf_parm+Pf(ix,ix)
        end do
     end if

  end if
  
end subroutine covariance_inflation

!-------------------------------------------------------------------------
! Make multiplicative parameter |
!-------------------------------------------------------------------------

subroutine make_inf(inf_parm)

  use parameter
  implicit none

  integer i,j
  integer iinf
  double precision,intent(out) :: inf_parm(ninf)

  if(mult_inf)then !Multiplicative inflation
     if(ninf == 1)then
        inf_parm(ninf)=1.05d0
     else
        do iinf=1,ninf
           inf_parm(iinf)=1.0d0+dble(iinf-1)*mult_amp
        end do
     end if
  elseif(add_inf)then !Additive inflation
     do iinf=1,ninf
        !Direct
        inf_parm(iinf)=dble(iinf)*add_amp
        !Logalithm
        !        inf_parm(iinf)=add_amp*10.d0**(0.1d0*dble(iinf-1))
        
     end do
  else if(rtpp_inf)then !RTPP
     do iinf=1,ninf
        inf_parm(iinf)=rtpp_intv*dble(iinf-1)
     end do
  else if(rtps_inf)then !RTPS
     do iinf=1,ninf
        inf_parm(iinf)=rtps_intv*dble(iinf-1)
     end do
  else
     write(*,*) "***No inflation***"
  end if
  
end subroutine make_inf
