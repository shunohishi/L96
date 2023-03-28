subroutine obs_operator(idt,nobs,H,HT)

  use parameter
  implicit none

  integer iobs
  integer opos(nobs)                                     !Obs. position

  integer,intent(in) :: idt ,nobs
  double precision,intent(out) :: H(nobs,nx),HT(nx,nobs) !Obs. operator

  if(nx == nobs)then
     do iobs=1,nobs
        opos(iobs)=iobs
     end do
  else
     !---Setting observation position
     call make_int_random(idt,1,nx,nobs,opos)
     call int_sort(nobs,opos)
  end if
     
  !---Observation Operator
  H(:,:)=0.d0
  do iobs=1,nobs
     H(iobs,opos(iobs))=1.d0
  end do
  HT(:,:)=transpose(H(:,:))
  
end subroutine obs_operator
