!----------------------------------------------------------
! Statistical value |
!-----------------------------------------------------------

subroutine ave_std(nx,x,ave,std)

  implicit none

  double precision,parameter :: dmiss=0.d0
  
  integer ix
  integer,intent(in) :: nx
  
  double precision,intent(in) :: x(nx)
  double precision,intent(out) :: ave,std
  double precision pass,miss
  
  !Average
  ave=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     ave=ave+x(ix)/dble(nx)
     pass=pass+1.d0
  end do
  
  !Standard deviation
  std=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     std=std+(x(ix)-ave)**2.d0/dble(nx-1)
     pass=pass+1.d0
  end do

  std=sqrt(std)
  
end subroutine ave_std

!----------------------------------

subroutine ave_std2(nx,ny,x,ave,std)

  implicit none

  integer ix,iy
  integer,intent(in) :: nx,ny
  
  double precision,intent(in) :: x(nx,ny)
  double precision,intent(out) :: ave,std
  
  !Average
  ave=0.d0
  do iy=1,ny
     do ix=1,nx
        ave=ave+x(ix,iy)
     end do
  end do
  ave=ave/dble(nx*ny)

  !Standard deviation
  std=0.d0
  do iy=1,ny
     do ix=1,nx
        std=std+(x(ix,iy)-ave)**2.d0
     end do
  end do
  std=sqrt(std/dble(nx*ny-1.d0))
  
end subroutine ave_std2

!--------------------------------------

subroutine ave_std_rmse(nx,x,ave,std)

  implicit none
  
  integer ix
  integer,intent(in) :: nx
  
  double precision,intent(in) :: x(nx) !spatially averaged RMSE 
  double precision,intent(out) :: ave,std
  double precision pass,miss
  
  !Average
  ave=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     ave=ave+x(ix)*x(ix)
     pass=pass+1.d0
  end do

  ave=sqrt(ave/dble(nx))
  
  !Standard deviation
  std=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     std=std+(x(ix)-ave)**2.d0/dble(nx-1)
     pass=pass+1.d0
  end do

  std=sqrt(std)
  
end subroutine ave_std_rmse

!--------------------------------------

subroutine bias_rmse(nx,x1,x2,bias,rmse)

  implicit none

  integer ix
  integer,intent(in) :: nx
  
  double precision,intent(in) :: x1(nx),x2(nx)
  double precision,intent(out) :: bias,rmse

  !Bias
  bias=0.d0
  do ix=1,nx
     bias=bias+(x1(ix)-x2(ix))/dble(nx)
  end do

  !RMSE
  rmse=0.d0
  do ix=1,nx
     rmse=rmse+(x1(ix)-x2(ix))**2.d0/dble(nx)
  end do
  rmse=sqrt(rmse)

end subroutine bias_rmse

!-----------------------------------------

subroutine correlation(n,x,y,cor)

  implicit none

  integer i
  integer,intent(in) :: n
  double precision,intent(in) :: x(n),y(n)
  double precision xmean,ymean
  double precision xsd,ysd
  double precision xycov
  double precision,intent(out) :: cor

  !Mean
  xmean=0.d0
  ymean=0.d0

  do i=1,n
     xmean=xmean+x(i)
     ymean=ymean+y(i)
  end do

  xmean=xmean/dble(n)
  ymean=ymean/dble(n)

  !Starndard deviation/covariance
  xsd=0.d0
  ysd=0.d0
  xycov=0.d0

  do i=1,n
     xsd=xsd+(x(i)-xmean)*(x(i)-xmean)
     ysd=ysd+(y(i)-ymean)*(y(i)-ymean)
     xycov=xycov+(x(i)-xmean)*(y(i)-ymean)
  end do

  xsd=sqrt(xsd/dble(n))
  ysd=sqrt(ysd/dble(n))
  xycov=xycov/dble(n)

  cor=xycov/(xsd*ysd)

end subroutine correlation

!---------------------------------------------------------

subroutine standarization(nx,x)

  implicit none

  integer ix
  integer,intent(in) :: nx

  double precision,intent(inout) :: x(nx)
  double precision ave,std

  call ave_std(nx,x,ave,std)

  do ix=1,nx
     x(ix)=(x(ix)-ave)/std
  end do
  
end subroutine standarization

!----------------------------------------------------------

subroutine diag_ave_std(nx,A,ave)

  implicit none

  double precision,parameter :: dmiss=0.d0
  
  integer ix
  integer,intent(in) :: nx

  double precision,intent(in) :: A(nx,nx)
  double precision,intent(out) :: ave
  double precision std
  double precision pass,miss

  !Average
  ave=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     if(A(ix,ix) == dmiss)then
        miss=miss+1.d0
     else
        ave=ave+A(ix,ix)/dble(nx)
        pass=pass+1.d0
     end if
  end do

  if(miss > 0.d0) ave=dmiss

  !Standard deviation
  std=0.d0
  pass=0.d0
  miss=0.d0
  do ix=1,nx
     if(A(ix,ix) == dmiss)then
        miss=miss+1.d0
     else
        std=std+(A(ix,ix)-ave)**2.d0/dble(nx-1)
        pass=pass+1.d0
     end if
  end do

  if(miss > 0.d0)then
     std=dmiss
  else
     std=sqrt(std)
  end if
  
end subroutine diag_ave_std

!----------------------------------------------------------

subroutine diag_off_ave_std(nx,A,diag_ave,diag_std,off_ave,off_std)

  implicit none

  integer ix,jx
  
  integer,intent(in) :: nx
  double precision,intent(in) :: A(nx,nx)

  double precision,intent(out) :: diag_ave,diag_std !Diagonal element
  double precision,intent(out) :: off_ave,off_std   !Off-diagonal element

  diag_ave=0.d0
  off_ave=0.d0
  
  do jx=1,nx
     do ix=1,nx
        if(ix == jx)then
           diag_ave=diag_ave+A(ix,jx)
        else
           off_ave=off_ave+A(ix,jx)
        end if
     end do
  end do

  diag_ave=diag_ave/dble(nx)
  off_ave=off_ave/dble(nx*nx-nx)

  diag_std=0.d0
  off_std=0.d0

  do jx=1,nx
     do ix=1,nx
        if(ix == jx)then
           diag_std=diag_std+(A(ix,jx)-diag_ave)**2.d0
        else
           off_std=off_std+(A(ix,jx)-off_ave)**2.d0
        end if
     end do
  end do

  diag_std=sqrt(diag_std/dble(nx-1))
  off_std=sqrt(off_std/dble(nx*nx-nx-1))
  
end subroutine diag_off_ave_std
