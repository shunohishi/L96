!------------------------------------------------------------------
! Make ensemble size |
!---------------------
!
!
!
!------------------------------------------------------------------
subroutine make_ensemble_size(nens,size)

  implicit none
  
  integer iens
  
  integer,intent(in) :: nens
  integer,intent(out) :: size(nens)

  if(nens == 1)then !Test case
     size(nens)=40
!     size(nens)=512
!  else if(11 < nens)then 
!     do iens=1,10
!        size(iens)=100*iens
!     end do
!     do iens=11,nens
!        size(iens)=size(10)+1000*(iens-10)
!     end do
  else
     do iens=1,nens
!        size(iens)=2**(4+iens-1) !16,32,64,128,256,512,...
        size(iens)=10*iens   !10,20,30,...
     end do
  end if
  
end subroutine make_ensemble_size
