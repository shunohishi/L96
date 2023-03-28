subroutine Pf_localization(nx,sigma_loc,Pf)

  implicit none

  integer ix,jx

  integer,intent(in) :: nx
  integer,intent(in) :: sigma_loc
  
  double precision,intent(inout) :: Pf(nx,nx)
  double precision dist,dist_loc
  
  dist_loc=2.d0*sqrt(10.d0/3.d0)*sigma_loc !Gaspari and Cohn (1999)
  
  do jx=1,nx
     do ix=1,nx
        
        dist=min(abs(dble(ix-nx-jx)),abs(dble(ix-jx)))
        dist=min(abs(dble(ix+nx-jx)),dist)

!        if(dist > dist_loc)then
!           Pf(ix,jx) = 0.d0
!        else
        Pf(ix,jx)=Pf(ix,jx)*exp(-0.5d0*(dist/dble(sigma_loc))**2.d0)
!        end if
           
     end do
  end do
  
end subroutine Pf_localization

!---------------------------------------------------------------

subroutine K_localization(nx,sigma_loc,H,inc)

  implicit none

  integer ix
  integer iobs

  integer,intent(in) :: nx
  integer,intent(in) :: sigma_loc
  
  double precision,intent(in) :: H(nx)
  double precision,intent(inout) :: inc(nx)
  double precision dist,dist_loc

  dist_loc=2.d0*sqrt(10.d0/3.d0)*sigma_loc !Gaspari and Cohn (1999)
  
  do ix=1,nx
     if(H(ix) == 1.d0)then
        iobs=ix
        exit
     end if
  end do

  do ix=1,nx
     dist=min(abs(dble(ix-nx-iobs)),abs(dble(ix-iobs)))
     dist=min(abs(dble(ix+nx-iobs)),dist)
!     if(dist > dist_loc)then
!        inc(ix)=0.d0
!     else
     inc(ix)=inc(ix)*exp(-0.5d0*(dist/dble(sigma_loc))**2.d0)
!     end if
  end do
  
end subroutine K_localization
