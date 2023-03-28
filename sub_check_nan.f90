subroutine check_nan(xt,x,xf,flag)

  use parameter
  implicit none

  integer ix,flag
  double precision,intent(inout) :: xt(nx),x(nx),xf(nx)

  do ix=1,nx
     if(isnan(xt(ix)) .or. isnan(x(ix)) .or. isnan(xf(ix)) &
          &  .or. abs(xt(ix)) > 999.d0 .or. abs(x(ix)) > 999.d0 .or. abs(xf(ix)) > 999.d0)then
        xt(:)=0.
        x(:)=0.
        xf(:)=0.
        flag=-1
        exit
     end if
  end do

end subroutine check_nan

!-----------------------------------------------------------------

subroutine check_P(it,type,P,flag)

  use parameter
  implicit none

  integer ix
  integer miss

  integer,intent(in) :: it,type
  integer,intent(inout) :: flag
  double precision,intent(in) :: P(nx,nx)

  flag=0
  miss=0
  do ix=1,nx
     if(P(ix,ix) <= 0.d0)then
        flag=1
        miss=miss+1
        if(type == 1) write(*,*) "Negative Pf:",it,ix,P(ix,ix)
        if(type == 2) write(*,*) "Negative Pa:",it,ix,P(ix,ix)
     end if
  end do

  if(miss == 0)then
     return
  else
!     stop
  endif
  
end subroutine check_P
