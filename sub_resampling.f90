!-------------------------------------------------------------------------
! SU: Stochastic universal |
!-------------------------------------------------------------------------

subroutine stochastic_universal(it,nens,cw,E)

  use parameter

  implicit none

  integer iens,jens
  integer,intent(in) :: it,nens
  integer pnum(nens)                  !Particle number at ith ensemble member
  
  double precision,intent(in) :: cw(0:nens)    !Cumulative weight
  double precision random                      !Random number
  double precision u(nens)
  double precision,intent(out) :: E(nens,nens) !Transform matrix
  
  !Partition
  call make_original_random(it,1,random)
  do iens=1,nens
     u(iens)=(random+dble(iens)-1.d0)/dble(nens)
  end do
  
  !Particle number by deterministic algorithm
  pnum(:)=0
  do iens=1,nens     
     do jens=1,nens
        if(u(iens) <= cw(jens))then
           pnum(iens)=jens
           exit
        end if        
     end do     
  end do
  
  E(:,:)=0.d0
  do iens=1,nens
     E(pnum(iens),iens)=1.d0
  end do
  
end subroutine stochastic_universal

!----------------------------------------------------------------------
! MN: Multinominal |
!----------------------------------------------------------------------

subroutine multinominal(it,nens,cw,E)

  implicit none

  integer iens,jens
  integer,intent(in) :: it,nens
  integer pnum(nens)
  
  double precision,intent(in) :: cw(0:nens) !Cumulative weight
  double precision random
  double precision,intent(out) :: E(nens,nens)

  pnum(:)=0
  do iens=1,nens

     call make_original_random(it+iens,1,random)

     do jens=1,nens
        if(random <= cw(jens))then
           pnum(iens)=jens
           exit
        end if
     end do
     
  end do

  do iens=1,nens
     if(pnum(iens) < 1 .or. nens < pnum(iens))then
        write(*,*) iens,pnum(iens)
        stop
     end if
  end do
  
  E(:,:)=0.d0
  do iens=1,nens
     E(pnum(iens),iens)=1.d0
  end do
  
end subroutine multinominal

!--------------------------------------------------------------------
! Residual resampling |
!--------------------------------------------------------------------

subroutine residual_resampling(it,nens,w,E)

  implicit none

  integer i,iens,jens

  integer,intent(in) :: it,nens

  double precision,intent(in):: w(nens)        !Weight
  double precision,intent(out) :: E(nens,nens) !Transform  matrix

  integer n(nens) !Number of particle: integer(Nw)
  integer ires,res
  double precision ratio(nens)
  double precision cum_ratio(0:nens)
  double precision random

  !The number of particle
  do iens=1,nens
     n(iens)=int(nens*w(iens))
  end do

  !Residual
  res=nens-sum(n)

  !Adjust residual
  if(res /= 0)then
     cum_ratio(0)=0.d0
     do iens=1,nens
        ratio(iens)=(dble(nens)*w(iens)-dble(n(iens)))/dble(res)
        cum_ratio(iens)=cum_ratio(iens-1)+ratio(iens)
     end do
     do ires=1,res
        call make_original_random(it+ires,1,random)
        do iens=1,nens
           if(cum_ratio(iens-1) <= random .and. random < cum_ratio(iens))then
              n(iens)=n(iens)+1
              exit
           end if
        end do
     end do
  end if

  !Make E matrix
  jens=0
  E(:,:)=0.d0
  do iens=1,nens
     if(n(iens) == 0) cycle
     do i=1,n(iens)
        jens=jens+1
        E(iens,jens)=1.d0
     end do
  end do
  
end subroutine residual_resampling
